use crate::database_setup::get_connection;
use super::{table::check_if_table_exists, Args};
use std::collections::HashMap;
use maplit::hashmap;

struct Adduct {
    adduct: String
}

fn get_adducts(matrix: String, prefix: &String) -> Vec<String> {
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

pub fn build_query(args: &Args, min_mz: f64, max_mz: f64, count: bool, prefix: String) -> String {
    check_if_table_exists("adducts", "user_adducts").unwrap();  
    let tissue_map: HashMap<&str, &str> = hashmap!{
        "HMDB (CSF)" => "csf",
        "HMDB (Urine)" => "urine",
        "HMDB (Serum)" => "serum",
    };

    let conventional_matrix: bool = args.matrix == "Positive Mode".to_string() || args.matrix == "Negative Mode".to_string();
    let endo_table = format!("{prefix}endogeneity", prefix=prefix);

    let met_type_string: String = build_query_with_table(&args.met_type,  &endo_table).unwrap();

    let mut query: String = String::new();

    query += "WITH concat_adducts AS (SELECT * FROM adducts UNION ALL SELECT * FROM user_adducts) ";
    let fg_or_adduct: String = get_adducts(args.matrix.clone(), &prefix).join(", ");

    if !count {
        if args.metabolome.starts_with("HMDB") {
            
            if !conventional_matrix {
                query += &format!(r#"SELECT DISTINCT (CAST({prefix}metabolites.mz AS REAL) + CASE WHEN concat_adducts.adduct IN ({fg_or_adduct}) THEN CAST(concat_adducts.deltamass AS REAL) ELSE 0 END) AS adjusted_mz, "#, fg_or_adduct=fg_or_adduct, prefix=prefix);
                query += &format!("{prefix}metabolites.name, concat_adducts.adduct, {prefix}db_accessions.hmdb, {prefix}metabolites.smiles, {prefix}metabolites.chemicalformula, ", prefix=prefix);
                query += &coverage_string(&build_condition_query3(&args.adducts, &prefix).unwrap(), &args.matrix, &prefix);
                query += " FROM";
                println!("cov_string: {:?}", &coverage_string(&build_condition_query3(&args.adducts, &prefix).unwrap(), &args.matrix, &prefix));
            } else {
                query += &format!(r#"SELECT (CAST({prefix}metabolites.mz AS REAL) + CASE WHEN concat_adducts.adduct IN ({fg_or_adduct}) THEN CAST(concat_adducts.deltamass AS REAL) ELSE 0 END) AS adjusted_mz, "#, fg_or_adduct=parse_fgs(&args.adducts), prefix=prefix);
                let num_adducts = args.adducts.len();
                query += &format!("{prefix}metabolites.name, concat_adducts.adduct, {prefix}db_accessions.hmdb, {prefix}metabolites.smiles, {prefix}metabolites.chemicalformula, {num_adducts} FROM ", num_adducts=num_adducts, prefix=prefix);
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
            query += &format!(" {prefix}metabolites ", prefix=prefix);
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
        cross_join_args = build_select_query(&get_adducts(args.matrix.clone(), &prefix)).unwrap_or("".to_string());
    } else {
        cross_join_args = build_select_query(&get_fgs(&args.adducts)).unwrap_or("".to_string());
    };
    if args.metabolome.starts_with("HMDB") {
        if !conventional_matrix {
            query += &format!("INNER JOIN {prefix}db_accessions ON {prefix}metabolites.id = {prefix}db_accessions.id INNER JOIN {prefix}endogeneity ON {prefix}metabolites.id = {prefix}endogeneity.id ");
            query += &format!("INNER JOIN {prefix}functional_groups ON {prefix}metabolites.id = {prefix}functional_groups.id INNER JOIN {prefix}derivatized_by ON {prefix}metabolites.id = {prefix}derivatized_by.id ");
            query += &format!("INNER JOIN {prefix}in_tissue ON {prefix}metabolites.id = {prefix}in_tissue.id ");
        } else {
            query += &format!("INNER JOIN {prefix}db_accessions ON {prefix}metabolites.id = {prefix}db_accessions.id INNER JOIN {prefix}endogeneity ON {prefix}metabolites.id = {prefix}endogeneity.id INNER JOIN {prefix}in_tissue ON {prefix}metabolites.id = {prefix}in_tissue.id ");
        }
    } 
    query += &format!(r#"CROSS JOIN ({a}) AS m LEFT JOIN concat_adducts ON m.mname = concat_adducts.adduct LEFT JOIN user_adducts ON m.mname = concat_adducts.adduct "#, a=&cross_join_args);


    

    let a: String = if !conventional_matrix {
        format!("WHERE concat_adducts.numfunctionalgroups <= {} AND adjusted_mz < {} AND adjusted_mz > {} {} AND ({})", &build_condition_query3(&args.adducts, &prefix).unwrap(), &max_mz.to_string(), &min_mz.to_string(), db_string, met_type_string)
    } else {
        if args.metabolome.starts_with("HMDB") {
            format!(r#"WHERE ({met_type}) AND CAST({prefix}metabolites.mz AS REAL) + CAST(concat_adducts.deltamass AS REAL) < {max} AND CAST({prefix}metabolites.mz AS REAL) + CAST(concat_adducts.deltamass AS REAL) > {min} {tissue}"#, 
            max=&max_mz.to_string(), min=&min_mz.to_string(), met_type=met_type_string, tissue=tissue_string, prefix=prefix)
        } else {
            format!(r#"WHERE CAST(lipids.mz AS REAL) + CAST(concat_adducts.deltamass AS REAL) < {max} AND CAST(lipids.mz AS REAL) + CAST(concat_adducts.deltamass AS REAL) > {min}"#, 
            max=&max_mz.to_string(), min=&min_mz.to_string())
        }

    };
    query += &a;
    if !count {
        if args.metabolome.starts_with("HMDB") {
            query += &format!(" ORDER BY CAST({prefix}metabolites.mz AS REAL) + CAST(concat_adducts.deltamass AS REAL)", prefix=prefix);
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

pub fn build_condition_query(types: &[String], start_with_and: bool, prefix: String) -> Option<String> {
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
        .map(|field| format!("{prefix}functional_groups.'{field}'", field=field, prefix=prefix)))
        .collect::<Vec<String>>()
        .join(" + ");

    query += " > 0)";
    Some(query)
}

pub fn build_condition_query3(types: &[String], prefix: &String) -> Option<String> {
    if types.is_empty() {
        return Some("".into());
    }

    let mut query: String ="(".to_string();

    query += &types.iter()
        .filter_map(|t| Some(t)
        .map(|field| format!("{prefix}functional_groups.'{field}'", field=field, prefix=prefix)))
        .collect::<Vec<String>>()
        .join(" + ");

    query += ")";
    Some(query)
}

pub fn build_condition_query2(types: &[String], start_with_and: bool, prefix: &String) -> Option<String> {
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
        .map(|field| format!("{prefix}endogeneity.'{field}' >= 1", field=field.to_lowercase(), prefix=prefix)))
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

pub fn coverage_string(fg_string: &str, matrix: &str, prefix: &String) -> String {
    let conn = get_connection().unwrap();
    let sql: &str = &format!("SELECT numfunctionalgroups, maxcoverage FROM {}adducts WHERE mname = '{}'", prefix, matrix);
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
        println!("adduct: {:?}", adduct);
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

    println!("prefix: {:?}\nfg_vec: {:?}", prefix, fg_vec);
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
    
    let mut query: String = "SELECT COUNT(*) FROM ".to_string();
    if met.starts_with("HMDB") {
        query += "metabolites INNER JOIN endogeneity ON metabolites.id = endogeneity.id INNER JOIN functional_groups ON metabolites.id = functional_groups.id INNER JOIN in_tissue ON metabolites.id = in_tissue.id";
    } else {
        query += "lipids";
        return query;
    };
    if matrix == "Positive Mode".to_string() || matrix == "Negative Mode".to_string() {
        query += &format!(" WHERE {met_type} {tissue_string}", 
            met_type= build_condition_query2(&typ, false, &"".to_string()).unwrap(),
            tissue_string = tissue_string)[..];
        return query;
    } else {
        query += &format!(" WHERE {met_type} {functional_group} {tissue_string}", 
            met_type = build_condition_query2(&typ, false, &"".to_string()).unwrap(),
            functional_group = build_condition_query(&adducts, true, "".to_string()).unwrap(),
            tissue_string = tissue_string)[..];
        return query;
    }
}

