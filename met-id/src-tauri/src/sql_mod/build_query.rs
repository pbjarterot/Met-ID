use crate::database_setup::get_connection;
use super::{table::check_if_table_exists, Args};
use std::collections::HashMap;
use maplit::hashmap;

struct Adduct {
    adduct: String
}

fn get_adducts(matrix: String) -> Vec<String> {
    let conn = get_connection().unwrap();
    let adduct_query = format!("WITH concat_adducts AS (SELECT * FROM adducts UNION ALL SELECT * FROM user_adducts) SELECT adduct FROM concat_adducts WHERE mname IN ('{}')", matrix);
    let mut stmt = conn.prepare(&adduct_query[..]).expect("Query cannot be run");

    let mut adducts: Vec<String> = Vec::new();
    let db_iter = stmt.query_map([], |row: &rusqlite::Row<'_>| {
        Ok(Adduct {
            adduct: row.get(0).unwrap_or("".to_string()),
        })
    }).unwrap();

    for (index, item) in db_iter.enumerate() {
        let row: Adduct = item.unwrap();
        adducts.insert(index, format!("'{}'", row.adduct));
    }

    adducts
}

fn get_fgs(fgs: &Vec<String>) -> Vec<String> {
    let mut new_vec = Vec::new();
    for index in 0..fgs.len() {
        new_vec.insert(index, format!("'{}'", fgs[index]))
    }
    new_vec
}

fn parse_fgs(fgs: &Vec<String>) -> String {
    let mut res_string = "".to_string();
    for i in 0..fgs.len() - 1 {
        res_string += &format!("'{}', ", fgs[i]);
    }
    res_string += &format!("'{}'", fgs[fgs.len()-1]);

    res_string
}

pub fn build_query(args: &Args, min_mz: f64, max_mz: f64, count: bool) -> String {
    check_if_table_exists("adducts", "user_adducts").unwrap();  
    let tissue_map: HashMap<&str, &str> = hashmap!{
        "HMDB (CSF)" => "csf",
        "HMDB (Urine)" => "urine",
        "HMDB (Serum)" => "serum",
    };

    

    let mut query: String = String::new();

    query += "WITH concat_adducts AS (SELECT * FROM adducts UNION ALL SELECT * FROM user_adducts), ";
    query += "concat_db_accessions AS (SELECT * FROM db_accessions UNION ALL SELECT * FROM user_db_accessions), ";
    query += "concat_metabolites AS (SELECT * FROM metabolites UNION ALL SELECT * FROM user_metabolites), ";
    query += "concat_derivatized_by AS (SELECT * FROM derivatized_by UNION ALL SELECT * FROM user_derivatized_by), ";
    query += "concat_endogeneity AS (SELECT * FROM endogeneity UNION ALL SELECT * FROM user_endogeneity), ";
    query += "concat_functional_groups AS (SELECT * FROM functional_groups UNION ALL SELECT * FROM user_functional_groups), ";
    query += "concat_in_tissue AS (SELECT * FROM in_tissue UNION ALL SELECT * FROM user_in_tissue) ";
    let fg_or_adduct: String = get_adducts(args.matrix.clone()).join(", ");

    let conventional_matrix: bool = args.matrix == "Positive Mode".to_string() || args.matrix == "Negative Mode".to_string();
    let endo_table: String = format!("concat_endogeneity");

    let met_type_string: String = build_query_with_table(&args.met_type,  &endo_table).unwrap();

    if !count {
        if args.metabolome.starts_with("HMDB") {
            
            if !conventional_matrix {
                query += &format!(r#"SELECT DISTINCT (CAST(concat_metabolites.mz AS REAL) + CASE WHEN concat_adducts.adduct IN ({fg_or_adduct}) THEN CAST(concat_adducts.deltamass AS REAL) ELSE 0 END) AS adjusted_mz, "#, fg_or_adduct=fg_or_adduct);
                query += &format!("concat_metabolites.name, concat_adducts.adduct, concat_db_accessions.hmdb, concat_metabolites.smiles, concat_metabolites.chemicalformula, ");
                query += &coverage_string(&build_condition_query3(&args.adducts).unwrap(), &args.matrix);
                query += " FROM";
            } else {
                query += &format!(r#"SELECT (CAST(concat_metabolites.mz AS REAL) + CASE WHEN concat_adducts.adduct IN ({fg_or_adduct}) THEN CAST(concat_adducts.deltamass AS REAL) ELSE 0 END) AS adjusted_mz, "#, fg_or_adduct=parse_fgs(&args.adducts));
                let num_adducts = args.adducts.len();
                query += &format!("concat_metabolites.name, concat_adducts.adduct, concat_db_accessions.hmdb, concat_metabolites.smiles, concat_metabolites.chemicalformula, {num_adducts} FROM ", num_adducts=num_adducts);
            }
        } else {
            query += &format!(r#"SELECT (CAST(lipids.mz AS REAL) + CASE WHEN concat_adducts.adduct IN ({fg_or_adduct}) THEN CAST(concat_adducts.deltamass AS REAL) ELSE 0 END) AS adjusted_mz, "#, fg_or_adduct=fg_or_adduct);
            query += &format!("lipids.name, concat_adducts.adduct, lipids.smiles, lipids.smiles, lipids.formula FROM");
        }
        
    } else {
        query += "SELECT COUNT(*) FROM";
    }

    let db_string: String;

    match args.metabolome.starts_with("HMDB") {
        true => {
            query += " concat_metabolites ";
            match tissue_map.get(&args.metabolome[..]){
                Some(x) => {
                    db_string = format!("AND (in_tissue.'{}' > 0)", x);
                },
                None => db_string = "".to_string(),
            }

        },
        false => {
            query += " lipids ";
            db_string = String::from("");
        },
    }
    let tissue_string:String = match tissue_map.get(&args.metabolome[..]){
        Some(x) => {
            format!("AND (in_tissue.'{}' >0 )", x)
        },
        None => "".to_string(),
    };
    
    let cross_join_args;
    if !conventional_matrix {
        cross_join_args = build_select_query(&get_adducts(args.matrix.clone())).unwrap_or("".to_string());
    } else {
        cross_join_args = build_select_query(&get_fgs(&args.adducts)).unwrap_or("".to_string());
    };
    if args.metabolome.starts_with("HMDB") {
        if !conventional_matrix {
            query += &format!("INNER JOIN concat_db_accessions ON concat_metabolites.id = concat_db_accessions.id INNER JOIN concat_endogeneity ON concat_metabolites.id = concat_endogeneity.id ");
            query += &format!("INNER JOIN concat_functional_groups ON concat_metabolites.id = concat_functional_groups.id INNER JOIN concat_derivatized_by ON concat_metabolites.id = concat_derivatized_by.id ");
            query += &format!("INNER JOIN concat_in_tissue ON concat_metabolites.id = concat_in_tissue.id ");
        } else {
            query += &format!("INNER JOIN concat_db_accessions ON concat_metabolites.id = concat_db_accessions.id INNER JOIN concat_endogeneity ON concat_metabolites.id = concat_endogeneity.id INNER JOIN concat_in_tissue ON concat_metabolites.id = concat_in_tissue.id ");
        }
    } 
    query += &format!(r#"CROSS JOIN ({a}) AS m LEFT JOIN concat_adducts ON m.mname = concat_adducts.adduct LEFT JOIN user_adducts ON m.mname = concat_adducts.adduct "#, a=&cross_join_args);


    

    let a: String = if !conventional_matrix {
        format!("WHERE concat_adducts.numfunctionalgroups <= {} AND adjusted_mz < {} AND adjusted_mz > {} {} AND ({})", &build_condition_query3(&args.adducts).unwrap(), &max_mz.to_string(), &min_mz.to_string(), db_string, met_type_string)
    } else {
        if args.metabolome.starts_with("HMDB") {
            format!(r#"WHERE ({met_type}) AND CAST(concat_metabolites.mz AS REAL) + CAST(concat_adducts.deltamass AS REAL) < {max} AND CAST(concat_metabolites.mz AS REAL) + CAST(concat_adducts.deltamass AS REAL) > {min} {tissue}"#, 
            max=&max_mz.to_string(), min=&min_mz.to_string(), met_type=met_type_string, tissue=tissue_string)
        } else {
            format!(r#"WHERE CAST(lipids.mz AS REAL) + CAST(concat_adducts.deltamass AS REAL) < {max} AND CAST(lipids.mz AS REAL) + CAST(concat_adducts.deltamass AS REAL) > {min}"#, 
            max=&max_mz.to_string(), min=&min_mz.to_string())
        }

    };
    query += &a;
    if !count {
        if args.metabolome.starts_with("HMDB") {
            query += &format!(" ORDER BY CAST(concat_metabolites.mz AS REAL) + CAST(concat_adducts.deltamass AS REAL)");
        } else {
            query += &format!(" ORDER BY CAST(lipids.mz AS REAL) + CAST(concat_adducts.deltamass AS REAL)");
        }
    }
    query
}

pub fn build_query_with_table(types: &[String], table: &str) -> Option<String> {
    if types.is_empty() {
        return Some("".to_string());
    }

    let query = types.iter()
        .filter_map(|t| Some(t).map(|field| format!("{}.'{}' >= 1", table, field.to_lowercase())))
        .collect::<Vec<String>>()
        .join(" OR ");

    Some(query)
}

pub fn build_condition_query(types: &[String], start_with_and: bool) -> Option<String> {
    if types.is_empty() {
        return Some("".into());
    }

    let mut query = if start_with_and {
        " AND (".to_string()
    } else {
        "(".to_string()
    };

    query += &types.iter()
        .filter_map(|t| Some(t)
        .map(|field| format!("concat_functional_groups.'{field}'", field=field)))
        .collect::<Vec<String>>()
        .join(" + ");

    query += " > 0)";
    Some(query)
}

pub fn build_condition_query3(types: &[String]) -> Option<String> {
    if types.is_empty() {
        return Some("".into());
    }

    let mut query: String ="(".to_string();

    query += &types.iter()
        .filter_map(|t| Some(t)
        .map(|field| format!("concat_functional_groups.'{field}'", field=field)))
        .collect::<Vec<String>>()
        .join(" + ");

    query += ")";
    Some(query)
}

pub fn build_condition_query2(types: &[String], start_with_and: bool) -> Option<String> {
    if types.is_empty() {
        return Some("".into());
    }

    let mut query = if start_with_and {
        " AND (".to_string()
    } else {
        "(".to_string()
    };

    query += &types.iter()
        .filter_map(|t| Some(t)
        .map(|field| format!("concat_endogeneity.'{field}' >= 1", field=field.to_lowercase())))
        .collect::<Vec<String>>()
        .join(" OR ");

    query += ")";
    Some(query)
}

pub fn build_select_query(types: &[String]) -> Option<String> {
    if types.is_empty() {
        return None;
    }
    
    let query = types.iter().enumerate()
        .filter_map(|(index, t)| Some(t).map(|field| {
            if index == 0 {
                format!("SELECT {} AS mname", field)
            } else {
                format!("UNION SELECT {}", field)
            }
        }))
        .collect::<Vec<String>>()
        .join(" ");

    Some(query)
}

pub fn coverage_string(fg_string: &str, matrix: &str) -> String {
    let conn = get_connection().unwrap();
    let sql: &str = &format!("WITH concat_adducts AS (SELECT * FROM adducts UNION ALL SELECT * FROM user_adducts) SELECT numfunctionalgroups, maxcoverage FROM concat_adducts WHERE mname = '{}'", matrix);
    let mut stmt = conn.prepare(sql).unwrap();
    let adducts_iter = stmt.query_map([], |row| {
        Ok((
            row.get::<_,i32>(0).unwrap(),
        row.get::<_,i32>(1).unwrap(),
    ))
    }).unwrap();

    let mut fg_vec: Vec<i32> = Vec::new();
    let mut cov_vec: Vec<i32> = Vec::new();

    for adduct in adducts_iter {
        match adduct {
            Ok((fg, cov)) => {
            if !fg_vec.contains(&fg) {
                fg_vec.push(fg);
                cov_vec.push(cov);
            }
            },
            Err(e) => println!("Error: {}", e)
        }
    }

    let max_index: usize = match fg_vec.iter()
                        .enumerate()
                        .max_by_key(|&(_, value)| value)
                        .map(|(index, _)| index) {
                            Some(val) => val,
                            None => 0
                        };
                        
    if max_index == 0 {
        return "0".to_string()
    } 

    let mut sql = "CASE ".to_string();

    sql += &format!("WHEN {} > {} THEN {} ", fg_string, fg_vec[max_index], cov_vec[max_index]);

    for index in 0..fg_vec.len() {
    sql += &format!("WHEN {} = {} THEN {} ", fg_string, fg_vec[index], cov_vec[index])
    }
    sql += "END as coverage_value";
    sql
}

pub fn build_count_query(met: String, matrix: String, typ: Vec<String>, adducts: Vec<String>) -> String {
    let tissue_map: HashMap<&str, &str> = hashmap!{
        "HMDB (CSF)" => "csf",
        "HMDB (Urine)" => "urine",
        "HMDB (Serum)" => "serum",
    };
    
    let tissue_string: String = match tissue_map.get(&met[..]){
        Some(x) => {
            format!("AND (in_tissue.'{}' > 0)", x)
        },
        None => "".to_string(),
    };
    
    let mut query: String = "WITH concat_endogeneity AS (SELECT * FROM endogeneity UNION ALL SELECT * FROM user_endogeneity), concat_functional_groups AS (SELECT * FROM functional_groups UNION ALL SELECT * FROM user_functional_groups) SELECT COUNT(*) FROM ".to_string();
    if met.starts_with("HMDB") {
        query += "metabolites INNER JOIN concat_endogeneity ON metabolites.id = concat_endogeneity.id INNER JOIN concat_functional_groups ON metabolites.id = concat_functional_groups.id INNER JOIN in_tissue ON metabolites.id = in_tissue.id";
    } else {
        query += "lipids";
        return query;
    };
    if matrix == "Positive Mode".to_string() || matrix == "Negative Mode".to_string() {
        query += &format!(" WHERE {met_type} {tissue_string}", 
            met_type= build_condition_query2(&typ, false).unwrap(),
            tissue_string = tissue_string)[..];
        return query;
    } else {
        query += &format!(" WHERE {met_type} {functional_group} {tissue_string}", 
            met_type = build_condition_query2(&typ, false).unwrap(),
            functional_group = build_condition_query(&adducts, true).unwrap(),
            tissue_string = tissue_string)[..];
        return query;
    }
}

