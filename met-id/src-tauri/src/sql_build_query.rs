use crate::sql::Args;
use std::collections::HashMap;


pub fn build_query(args: Args, min_mz: f64, max_mz: f64, count: bool) -> String {
    let mut met_type_map: HashMap<&str, &str> = HashMap::new();
    met_type_map.insert("Endogenous", "endogenous");
    met_type_map.insert("Exogenous", "exogenous");
    met_type_map.insert("Unspecified", "unspecified");

    /* 
    let _matrix_map: HashMap<&str, &str> = HashMap::new();
        "FMP-10" => "fmp",
        "DPP" => "dpp",
        "TAHS" => "tahs",
        "CA" => "ca",
        "Girard P" => "girardp",
        "Girard T" => "girardt",
        "AMPP" => "ampp",
        "Boronic Acid" => "boronicacid"
    */

    let mut functional_groups_map: HashMap<&str, &str> = HashMap::new();
    functional_groups_map.insert("Phenols", "phenols");
    functional_groups_map.insert("Catechols", "catechols");
    functional_groups_map.insert("Carbonyls", "carbonyls");
    functional_groups_map.insert("Aldehydes", "aldehydes");
    functional_groups_map.insert("Carboxylic", "carboxylicacids");
    functional_groups_map.insert("Primary", "primaryamines");
    functional_groups_map.insert("Triols", "triols");
    functional_groups_map.insert("Diols", "diols");
    functional_groups_map.insert("Hydroxyls", "hydroxyls");

    let mut adduct_map: HashMap<&str, &str> = HashMap::new();
    adduct_map.insert("[M+K]+","M+K");
    adduct_map.insert("[M+H]+" , "M+H");
    adduct_map.insert("[M+Na]+","M+Na");
    adduct_map.insert("[M+K-2H]-" ,"M+K-H2");
    adduct_map.insert("[M-H]-" , "M-H");
    adduct_map.insert("[M+Na-2H]-" , "M+Na-H2");
    adduct_map.insert("[M+Cl]-" ,"M+Cl");
    adduct_map.insert("M+FMP10" , "M+FMP10");
    adduct_map.insert("M+2FMP10a", "M+2FMP10a");
    adduct_map.insert("M+2FMP10b", "M+2FMP10b");
    adduct_map.insert("M+3FMP10a", "M+3FMP10a");
    adduct_map.insert("M+3FMP10b", "M+3FMP10b");
    adduct_map.insert("M+3FMP10c", "M+3FMP10c");
    //extend for more matrices


    let mut tissue_map: HashMap<&str, &str> = HashMap::new();
    tissue_map.insert("HMDB (CSF)", "csf");
    tissue_map.insert("HMDB (Urine)", "urine");
    tissue_map.insert("HMDB (Serum)", "serum");



    let if_fmp10: bool = args.matrix == "FMP-10".to_string();
    let pos_neg: bool = args.matrix == "Positive mode".to_string() || args.matrix == "Negative mode".to_string();
    let met_type_string: String = build_query_with_table(&args.met_type, &met_type_map, "endogeneity").unwrap();
    let chosen_deriv_string: &str = &args.matrix;

    let mut query: String = String::new();

    let fg_or_adduct: String = {
        if if_fmp10 {
            let adduct_args: Vec<String> = vec![String::from("M+FMP10"), String::from("M+2FMP10a"), String::from("M+2FMP10b"), String::from("M+3FMP10a"), String::from("M+3FMP10b"), String::from("M+3FMP10c")];
            build_custom_query(&adduct_args, &adduct_map, "adduct").unwrap_or("".to_string())
        } else if pos_neg {
            build_custom_query(&args.adducts, &adduct_map, "adduct").unwrap_or("".to_string())
        } else {
            build_custom_query(&args.adducts, &functional_groups_map, "fg").unwrap_or("".to_string())
        }
    };
    println!("FG OR ADDUCT: {:}", fg_or_adduct);


    if !count {
        //query += "SELECT (CAST(metabolites.mz AS REAL) + CAST(matrices.deltamass AS REAL)), ";
        query += &format!(r#"SELECT (CAST(metabolites.mz AS REAL) + CASE WHEN adducts.adduct IN ({fg_or_adduct}) THEN CAST(adducts.deltamass AS REAL) ELSE 0 END) AS adjusted_mz,"#, fg_or_adduct=fg_or_adduct);
        query += "metabolites.name, adducts.adduct, db_accessions.hmdb, metabolites.smiles, metabolites.chemicalformula FROM";
    } else {
        query += "SELECT COUNT(*) FROM";
    }

    let mut db_string: String;

    match args.metabolome.starts_with("HMDB") {
        true => {
            query += " metabolites ";
            match tissue_map.get(&args.metabolome[..]){
                Some(x) => {
                    db_string = "AND (".to_string();
                    db_string += x;
                    db_string += &String::from(" >0 )");
                },
                None => db_string = "".to_string(),
            }

        },
        false => {
            query += " lipids ";
            db_string = String::from("");
        },
    }


    let cross_join_args = {
        if if_fmp10 {
            let adduct_args: Vec<String> = vec![String::from("M+FMP10"), String::from("M+2FMP10a"), String::from("M+2FMP10b"), String::from("M+3FMP10a"), String::from("M+3FMP10b"), String::from("M+3FMP10c")];
            build_select_query(&adduct_args, &adduct_map).unwrap_or("".to_string())
        } else if pos_neg {
            build_select_query(&args.adducts, &adduct_map).unwrap_or("".to_string())
        } else {
            //extend this for more matrices
            build_select_query(&args.adducts, &adduct_map).unwrap_or("".to_string())
        }
    };

    //query += "a INNER JOIN (matrices b, db_accessions d, endogeneity e, functional_groups f) ";
    query += "INNER JOIN db_accessions ON metabolites.id = db_accessions.id INNER JOIN endogeneity ON metabolites.id = endogeneity.id INNER JOIN functional_groups ON metabolites.id = functional_groups.id INNER JOIN derivatized_by ON metabolites.id = derivatized_by.id ";
    query += &format!(r#"CROSS JOIN ({a}) AS m LEFT JOIN adducts ON m.mname = adducts.adduct "#, a=&cross_join_args);
    

    let a: String = if if_fmp10 {
        format!(r#"WHERE CASE WHEN adducts.adduct='M+FMP10' AND derivatized_by.fmp = 1 AND CAST(metabolites.mz AS REAL) + CAST(adducts.deltamass AS REAL) < {max} AND CAST(metabolites.mz AS REAL) + CAST(adducts.deltamass AS REAL) > {min} {db} AND ({met_type}) AND ({func_group}) THEN 1 WHEN adducts.adduct IN ('M+Anth', 'M+2Anth', 'M+2Anth_2') AND derivatized_by.fmp = 2 AND CAST(metabolites.mz AS REAL) + CAST(adducts.deltamass AS REAL) < {max} AND CAST(metabolites.mz AS REAL) + CAST(adducts.deltamass AS REAL) > {min} {db} AND ({met_type}) AND ({func_group}) THEN 1 WHEN adducts.adduct IN ('M+Anth', 'M+2Anth', 'M+2Anth_2', 'M+3Anth1', 'M+3Anth2', 'M+3Anth3') AND derivatized_by.fmp>2 AND CAST(metabolites.mz AS REAL) + CAST(adducts.deltamass AS REAL) < {max} AND CAST(metabolites.mz AS REAL) + CAST(adducts.deltamass AS REAL) > {min} {db} AND ({met_type}) AND ({func_group}) THEN 1 ELSE 0 END=1"#, 
        max=&max_mz.to_string(), min=&min_mz.to_string(), db=db_string, met_type = met_type_string, func_group = build_query_with_table(&args.adducts, &functional_groups_map, "functional_groups").unwrap())
    } else if pos_neg {
        format!(r#"WHERE ({met_type}) AND CAST(metabolites.mz AS REAL) + CAST(adducts.deltamass AS REAL) < {max} AND CAST(metabolites.mz AS REAL) + CAST(adducts.deltamass AS REAL) > {min} AND adducts.adduct IN ({adduct})"#, 
        max=&max_mz.to_string(), min=&min_mz.to_string(), 
        //db=db_string, 
        met_type=met_type_string, 
        adduct=build_csv_query(&args.adducts, &adduct_map).unwrap().replace("\\", ""))

    } else {
        format!(r#"WHERE adducts.adduct ='M+{deriv}' AND CAST(metabolites.mz AS REAL) + CAST(adducts.deltamass AS REAL) < {max} AND CAST(metabolites.mz AS REAL) + CAST(adducts.deltamass AS REAL) > {min} {db} AND ({met_type}) AND ({func_group})"#, 
        max=&max_mz.to_string(), 
        min=&min_mz.to_string(),
        db=db_string, 
        met_type = met_type_string, 
        func_group=build_query_with_table(&args.adducts, &functional_groups_map, "functional_groups").unwrap(),
        deriv = chosen_deriv_string
        )
    };
    query += &a;
    if !count {
        query += " ORDER BY CAST(metabolites.mz AS REAL) + CAST(adducts.deltamass AS REAL)";
    }
    query
}


pub fn build_query_with_table(types: &[String], mapping: &HashMap<&str, &str>, table: &str) -> Option<String> {
    if types.is_empty() {
        return None;
    }

    let query = types.iter()
        .filter_map(|t| mapping.get(t.as_str()).map(|&field| format!("{}.{} = 1", table, field)))
        .collect::<Vec<String>>()
        .join(" OR ");

    Some(query)
}

pub fn build_csv_query(types: &[String], mapping: &HashMap<&str, &str>) -> Option<String> {
    if types.is_empty() {
        return None;
    }
    let query = types.iter()
        .filter_map(|t| mapping.get(t.as_str()))
        .map(|&field| format!("'{}'", field))
        .collect::<Vec<String>>()
        .join(", ");

    Some(query)
}

pub fn build_condition_query(types: &[String], mapping: &HashMap<&str, &str>, start_with_and: bool) -> Option<String> {
    if types.is_empty() {
        return Some("".into());
    }

    let mut query = if start_with_and {
        " AND (".to_string()
    } else {
        "(".to_string()
    };

    query += &types.iter()
        .filter_map(|t| mapping.get(t.as_str()))
        .map(|&field| format!("{} > 0", field))
        .collect::<Vec<String>>()
        .join(" OR ");

    query += ")";
    Some(query)
}

pub fn build_select_query(types: &[String], mapping: &HashMap<&str, &str>) -> Option<String> {
    if types.is_empty() {
        return None;
    }
    
    let query = types.iter().enumerate()
        .filter_map(|(index, t)| mapping.get(t.as_str()).map(|&field| {
            if index == 0 {
                format!("SELECT '{}' AS mname", field)
            } else {
                format!("UNION SELECT '{}'", field)
            }
        }))
        .collect::<Vec<String>>()
        .join(" ");

    Some(query)
}

pub fn build_custom_query(types: &[String], mapping: &HashMap<&str, &str>, fg_or_adduct: &str) -> Option<String> {
    if types.is_empty() {
        return None;
    }

    let query = types.iter()
        .filter_map(|t| mapping.get(t.as_str()))
        .map(|&field| {
            if fg_or_adduct == "adduct" {
                format!("'{}'", field)
            } else if fg_or_adduct == "fg" {
                format!("functional_groups.{}", field)
            } else {
                "".to_string()
            }
        })
        .collect::<Vec<String>>()
        .join(", ");

    Some(query)
}