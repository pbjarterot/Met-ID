use super::{table::check_if_table_exists, Args};
use crate::database_setup::get_connection;
use maplit::hashmap;
use rusqlite::Result;
use std::collections::HashMap;

/* 
This file contains functions to build the SQLite query for the mass matching of the database


Patrik Bjärterot - March 2025
*/


struct Adduct {
    adduct: String,
}

fn get_adducts(matrix: String) -> Vec<String> {
    let conn = get_connection().unwrap();
    let adduct_query = format!("WITH temp_concat_adducts AS (SELECT * FROM adducts UNION ALL SELECT * FROM user_adducts) SELECT adduct FROM temp_concat_adducts WHERE mname IN ('{}')", matrix);
    let mut stmt = conn
        .prepare(&adduct_query[..])
        .expect("Query cannot be run");

    let mut adducts: Vec<String> = Vec::new();
    let db_iter = stmt
        .query_map([], |row: &rusqlite::Row<'_>| {
            Ok(Adduct {
                adduct: row.get(0).unwrap_or("".to_string()),
            })
        })
        .unwrap();

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
    res_string += &format!("'{}'", fgs[fgs.len() - 1]);

    res_string
}

fn table_exists(tx: &rusqlite::Transaction<'_>) -> Result<usize> {
    tx.execute("SELECT COUNT(*) FROM temp_concat_adducts", [])
}

pub fn check_temp_tables() {
    let mut conn: r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager> =
        get_connection().unwrap();

    check_if_table_exists("metabolites", "user_metabolites").unwrap();
    check_if_table_exists("derivatized_by", "user_derivatized_by").unwrap();
    check_if_table_exists("endogeneity", "user_endogeneity").unwrap();
    check_if_table_exists("in_tissue", "user_in_tissue").unwrap();
    check_if_table_exists("db_accessions", "user_db_accessions").unwrap();
    check_if_table_exists("functional_groups", "user_functional_groups").unwrap();
    check_if_table_exists("adducts", "user_adducts").unwrap();

    let tx: rusqlite::Transaction<'_> = conn.transaction().unwrap();
    match table_exists(&tx) {
        Ok(_) => (),
        Err(rusqlite::Error::ExecuteReturnedResults) => (),
        Err(_) => {
            tx.execute(
                "CREATE TEMP TABLE temp_concat_adducts AS
                SELECT * FROM adducts
                UNION ALL
                SELECT * FROM user_adducts",
                [],
            )
            .unwrap();

            tx.execute(
                "CREATE TEMP TABLE temp_concat_db_accessions AS
                SELECT * FROM db_accessions
                UNION ALL
                SELECT * FROM user_db_accessions",
                [],
            )
            .unwrap();

            tx.execute(
                "CREATE TEMP TABLE temp_concat_metabolites AS
                SELECT * FROM metabolites
                UNION ALL
                SELECT * FROM user_metabolites",
                [],
            )
            .unwrap();

            tx.execute(
                "CREATE TEMP TABLE temp_concat_lipids AS
                SELECT * FROM lipids 
                UNION ALL 
                SELECT * FROM user_metabolites",
                [],
            )
            .unwrap();

            tx.execute(
                "CREATE TEMP TABLE temp_concat_derivatized_by AS
                SELECT * FROM derivatized_by
                UNION ALL
                SELECT * FROM user_derivatized_by",
                [],
            )
            .unwrap();

            tx.execute(
                "CREATE TEMP TABLE temp_concat_endogeneity AS
                SELECT * FROM endogeneity
                UNION ALL
                SELECT * FROM user_endogeneity",
                [],
            )
            .unwrap();

            tx.execute(
                "CREATE TEMP TABLE temp_concat_functional_groups AS
                SELECT * FROM functional_groups
                UNION ALL
                SELECT * FROM user_functional_groups",
                [],
            )
            .unwrap();
            tx.execute(
                "CREATE TEMP TABLE temp_concat_lipid_functional_groups AS
                SELECT * FROM lipids_functional_groups 
                UNION ALL
                SELECT * FROM user_functional_groups",
                [],
            )
            .unwrap();

            tx.execute(
                "CREATE TEMP TABLE temp_concat_in_tissue AS
                SELECT * FROM in_tissue
                UNION ALL
                SELECT * FROM user_in_tissue",
                [],
            )
            .unwrap();

            tx.execute(
                "CREATE INDEX idx_temp_concat_db_accessions_id ON temp_concat_db_accessions(id)",
                [],
            )
            .unwrap();

            tx.execute(
                "CREATE INDEX idx_temp_concat_metabolites_id ON temp_concat_metabolites(id)",
                [],
            )
            .unwrap();

            tx.execute(
                "CREATE INDEX idx_temp_concat_lipids_id ON temp_concat_lipids(id)",
                [],
            )
            .unwrap();

            tx.execute(
                "CREATE INDEX idx_temp_concat_derivatized_by_id ON temp_concat_derivatized_by(id)",
                [],
            )
            .unwrap();

            tx.execute(
                "CREATE INDEX idx_temp_concat_endogeneity_id ON temp_concat_endogeneity(id)",
                [],
            )
            .unwrap();
            tx.execute("CREATE INDEX idx_temp_concat_endogeneity_criteria ON temp_concat_endogeneity(endogenous, exogenous, unspecified)", []).unwrap();

            tx.execute("CREATE INDEX idx_temp_concat_functional_groups_id ON temp_concat_functional_groups(id)", []).unwrap();
            tx.execute("CREATE INDEX idx_temp_concat_lipid_functional_groups_id ON temp_concat_lipid_functional_groups(id)", []).unwrap();

            tx.execute("CREATE INDEX idx_temp_concat_functional_groups_criteria ON temp_concat_functional_groups('Phenolic Hydroxyls', 'Primary Amines')", []).unwrap();

            tx.execute(
                "CREATE INDEX idx_temp_concat_in_tissue_id ON temp_concat_in_tissue(id)",
                [],
            )
            .unwrap();

            tx.commit().unwrap();
        }
    }
}

pub fn build_query(args: &Args, min_mz: f64, max_mz: f64, count: bool) -> String {
    check_if_table_exists("adducts", "user_adducts").unwrap();
    let tissue_map: HashMap<&str, &str> = hashmap! {
        "HMDB (CSF)" => "csf",
        "HMDB (Urine)" => "urine",
        "HMDB (Serum)" => "serum",
    };

    check_temp_tables();

    let mut query: String = String::new();
    let fg_or_adduct: String = get_adducts(args.matrix.clone()).join(", ");

    let conventional_matrix: bool =
        args.matrix == "Positive Mode".to_string() || args.matrix == "Negative Mode".to_string();
    let endo_table: String = format!("temp_concat_endogeneity");

    let met_type_string: String = build_query_with_table(&args.met_type, &endo_table).unwrap();

    let tissue_string: String = match tissue_map.get(&args.metabolome[..]) {
        Some(x) => {
            format!("AND (in_tissue.'{}' >0 )", x)
        }
        None => "".to_string(),
    };

    let cross_join_args;
    if !conventional_matrix {
        cross_join_args =
            build_select_query(&get_adducts(args.matrix.clone())).unwrap_or("".to_string());
    } else {
        cross_join_args = build_select_query(&get_fgs(&args.adducts)).unwrap_or("".to_string());
    };

    println!("\n\n\n{:?}\n\n\n", cross_join_args);
    let mut db_string: String = "".to_string();

    if args.metabolome.starts_with("HMDB") {
        if let Some(x) = tissue_map.get(&args.metabolome[..]) {
            db_string = format!("AND (in_tissue.'{}' > 0)", x);
        }
    };

    match count {
        true => {
            query += "SELECT COUNT(*) FROM";
        }
        false => match (args.metabolome.starts_with("HMDB"), conventional_matrix) {
            (true, false) => {
                //HMDB & Derivatizing matrix
                query += &format!(
                    "SELECT DISTINCT (temp_concat_metabolites.mz + CASE WHEN temp_concat_adducts.adduct IN ({fg_or_adduct}) \
                    THEN temp_concat_adducts.deltamass ELSE 0 END) \
                    AS adjusted_mz, temp_concat_metabolites.name, temp_concat_adducts.adduct, temp_concat_db_accessions.hmdb, \
                    temp_concat_metabolites.smiles, temp_concat_metabolites.chemicalformula, {coverage_string} FROM temp_concat_metabolites ",
                    fg_or_adduct = fg_or_adduct,
                    coverage_string = &coverage_string(&build_condition_query3(&args.adducts, "").unwrap(), &args.matrix)
                );
            }
            (true, true) => {
                //HMDB & Pos neg mode
                query += &format!(
                    "SELECT DISTINCT (temp_concat_metabolites.mz + CASE WHEN temp_concat_adducts.adduct IN ({fg_or_adduct}) \
                    THEN temp_concat_adducts.deltamass ELSE 0 END) \
                    AS adjusted_mz, temp_concat_metabolites.name, temp_concat_adducts.adduct, temp_concat_db_accessions.hmdb, \
                    temp_concat_metabolites.smiles, temp_concat_metabolites.chemicalformula, {num_adducts} FROM temp_concat_metabolites ",
                    fg_or_adduct = parse_fgs(&args.adducts),
                    num_adducts = args.adducts.len()
                );
            }
            (false, false) => {
                //lipidmaps & Derivatizing matrix
                query += &format!(
                    "SELECT DISTINCT (temp_concat_lipids.mz + CASE WHEN temp_concat_adducts.adduct IN ({fg_or_adduct}) \
                    THEN temp_concat_adducts.deltamass ELSE 0 END) \
                    AS adjusted_mz, temp_concat_lipids.name, temp_concat_adducts.adduct, '0', \
                    temp_concat_lipids.smiles, temp_concat_lipids.formula, {num_adducts} FROM temp_concat_lipids ",
                    fg_or_adduct = fg_or_adduct,
                    num_adducts = &coverage_string(&build_condition_query3(&args.adducts, "_lipid").unwrap(), &args.matrix)
                    );
            }
            (false, true) => {
                //lipidmaps & pos neg mode
                query += &format!(
                    "SELECT DISTINCT (temp_concat_lipids.mz + CASE WHEN temp_concat_adducts.adduct IN ({fg_or_adduct}) \
                    THEN temp_concat_adducts.deltamass ELSE 0 END) \
                    AS adjusted_mz, temp_concat_lipids.name, temp_concat_adducts.adduct, '0', \
                    temp_concat_lipids.smiles, temp_concat_lipids.formula, {num_adducts} FROM temp_concat_lipids ",
                    fg_or_adduct = parse_fgs(&args.adducts),
                    num_adducts = args.adducts.len()
                );
            }
        },
    }

    match (args.metabolome.starts_with("HMDB"), conventional_matrix) {
        (true, false) => {
            query += &format!("INNER JOIN temp_concat_db_accessions ON temp_concat_metabolites.id = temp_concat_db_accessions.id \
                                INNER JOIN temp_concat_endogeneity ON temp_concat_metabolites.id = temp_concat_endogeneity.id \
                                INNER JOIN temp_concat_functional_groups ON temp_concat_metabolites.id = temp_concat_functional_groups.id \
                                INNER JOIN temp_concat_derivatized_by ON temp_concat_metabolites.id = temp_concat_derivatized_by.id \
                                INNER JOIN temp_concat_in_tissue ON temp_concat_metabolites.id = temp_concat_in_tissue.id \
                                CROSS JOIN ({a}) AS m LEFT JOIN temp_concat_adducts ON m.mname = temp_concat_adducts.adduct \
                                LEFT JOIN user_adducts ON m.mname = temp_concat_adducts.adduct \
                                WHERE temp_concat_adducts.numfunctionalgroups <= {b} AND adjusted_mz < {max} AND adjusted_mz > {min} {db_string}", 
                                a = &cross_join_args,
                                b = &build_condition_query3(&args.adducts, "").unwrap(), 
                                max = &max_mz.to_string(), 
                                min = &min_mz.to_string(), 
                                db_string = db_string
                                );
        }
        (true, true) => {
            query += &format!("INNER JOIN temp_concat_db_accessions ON temp_concat_metabolites.id = temp_concat_db_accessions.id \
                                INNER JOIN temp_concat_endogeneity ON temp_concat_metabolites.id = temp_concat_endogeneity.id \
                                INNER JOIN temp_concat_in_tissue ON temp_concat_metabolites.id = temp_concat_in_tissue.id \
                                CROSS JOIN ({a}) AS m LEFT JOIN temp_concat_adducts ON m.mname = temp_concat_adducts.adduct \
                                LEFT JOIN user_adducts ON m.mname = temp_concat_adducts.adduct \
                                AND ({met_type_string}) WHERE ({met_type}) AND temp_concat_metabolites.mz + temp_concat_adducts.deltamass < {max} AND temp_concat_metabolites.mz + temp_concat_adducts.deltamass > {min} {tissue}",
                                a = cross_join_args,
                                max = &max_mz.to_string(),
                                min = &min_mz.to_string(),
                                met_type = met_type_string,
                                tissue = tissue_string
                            );
        }
        (false, false) => {
            query += &format!("INNER JOIN temp_concat_lipid_functional_groups ON temp_concat_lipids.id = temp_concat_lipid_functional_groups.id \
                            CROSS JOIN ({a}) AS m LEFT JOIN temp_concat_adducts ON m.mname = temp_concat_adducts.adduct \
                            LEFT JOIN user_adducts ON m.mname = temp_concat_adducts.adduct \
                            WHERE temp_concat_lipids.mz + temp_concat_adducts.deltamass < {max} AND temp_concat_lipids.mz + temp_concat_adducts.deltamass > {min}",
                            a = cross_join_args,
                            max = &max_mz.to_string(),
                            min = &min_mz.to_string());
        }
        (false, true) => {
            query += &format!("CROSS JOIN ({a}) AS m LEFT JOIN temp_concat_adducts ON m.mname = temp_concat_adducts.adduct \
                                LEFT JOIN user_adducts ON m.mname = temp_concat_adducts.adduct \
                                WHERE temp_concat_lipids.mz + temp_concat_adducts.deltamass < {max} AND temp_concat_lipids.mz + temp_concat_adducts.deltamass > {min}", 
                                a = cross_join_args,
                                max = &max_mz.to_string(),
                                min = &min_mz.to_string());
        }
    }

    if !count {
        query += &format!(" ORDER BY adjusted_mz");
    }
    query.trim().to_string()
}

pub fn build_query_with_table(types: &[String], table: &str) -> Option<String> {
    if types.is_empty() {
        return Some("".to_string());
    }

    let query = types
        .iter()
        .filter_map(|t| Some(t).map(|field| format!("{}.'{}' >= 1", table, field.to_lowercase())))
        .collect::<Vec<String>>()
        .join(" OR ");

    Some(query)
}

pub fn build_condition_query(
    types: &[String],
    start_with_and: bool,
    metabolome: &str,
) -> Option<String> {
    if types.is_empty() {
        return Some("".into());
    }

    let mut query = if start_with_and {
        " AND (".to_string()
    } else {
        "(".to_string()
    };

    query += &types
        .iter()
        .filter_map(|t| {
            Some(t).map(|field| {
                format!(
                    "temp_concat{metabolome}_functional_groups.'{field}'",
                    metabolome = metabolome,
                    field = field
                )
            })
        })
        .collect::<Vec<String>>()
        .join(" + ");

    query += " > 0)";
    Some(query)
}

pub fn build_condition_query3(types: &[String], metabolome: &str) -> Option<String> {
    if types.is_empty() {
        return Some("".into());
    }

    let mut query: String = "(".to_string();

    query += &types
        .iter()
        .filter_map(|t| {
            Some(t).map(|field| {
                format!(
                    "temp_concat{metabolome}_functional_groups.'{field}'",
                    metabolome = metabolome,
                    field = field
                )
            })
        })
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

    query += &types
        .iter()
        .filter_map(|t| {
            Some(t).map(|field| {
                format!(
                    "temp_concat_endogeneity.'{field}' >= 1",
                    field = field.to_lowercase()
                )
            })
        })
        .collect::<Vec<String>>()
        .join(" OR ");

    query += ")";
    Some(query)
}

pub fn build_select_query(types: &[String]) -> Option<String> {
    if types.is_empty() {
        return None;
    }

    let query = types
        .iter()
        .enumerate()
        .filter_map(|(index, t)| {
            Some(t).map(|field| {
                if index == 0 {
                    format!("SELECT {} AS mname", field)
                } else {
                    format!("UNION SELECT {}", field)
                }
            })
        })
        .collect::<Vec<String>>()
        .join(" ");

    Some(query)
}

pub fn coverage_string(fg_string: &str, matrix: &str) -> String {
    let conn = get_connection().unwrap();
    let sql: &str = &format!("WITH temp_concat_adducts AS (SELECT * FROM adducts UNION ALL SELECT * FROM user_adducts) SELECT numfunctionalgroups, maxcoverage FROM temp_concat_adducts WHERE mname = '{}'", matrix);
    let mut stmt = conn.prepare(sql).unwrap();
    let adducts_iter = stmt
        .query_map([], |row| {
            Ok((row.get::<_, i32>(0).unwrap(), row.get::<_, i32>(1).unwrap()))
        })
        .unwrap();

    let mut fg_vec: Vec<i32> = Vec::new();
    let mut cov_vec: Vec<i32> = Vec::new();

    for adduct in adducts_iter {
        match adduct {
            Ok((fg, cov)) => {
                if !fg_vec.contains(&fg) {
                    fg_vec.push(fg);
                    cov_vec.push(cov);
                }
            }
            Err(e) => println!("Error with the code: {}", e),
        }
    }
    let max_index: usize = match fg_vec
        .iter()
        .enumerate()
        .max_by_key(|&(_, value)| value)
        .map(|(index, _)| index)
    {
        Some(val) => val,
        None => 0,
    };
    /*
    if max_index == 0 {
        return "0".to_string()
    }
    */
    let mut sql = "CASE ".to_string();

    sql += &format!(
        "WHEN {} > {} THEN {} ",
        fg_string, fg_vec[max_index], cov_vec[max_index]
    );

    for index in 0..fg_vec.len() {
        sql += &format!(
            "WHEN {} = {} THEN {} ",
            fg_string, fg_vec[index], cov_vec[index]
        )
    }
    sql += "END as coverage_value";
    sql
}

pub fn build_count_query(
    met: String,
    matrix: String,
    typ: Vec<String>,
    adducts: Vec<String>,
) -> String {
    let tissue_map: HashMap<&str, &str> = hashmap! {
        "HMDB (CSF)" => "csf",
        "HMDB (Urine)" => "urine",
        "HMDB (Serum)" => "serum",
    };

    let tissue_string: String = match tissue_map.get(&met[..]) {
        Some(x) => {
            format!("AND (in_tissue.'{}' > 0)", x)
        }
        None => "".to_string(),
    };

    let mut query: String = "WITH temp_concat_endogeneity AS (SELECT * FROM endogeneity UNION ALL SELECT * FROM user_endogeneity), temp_concat_functional_groups AS (SELECT * FROM functional_groups UNION ALL SELECT * FROM user_functional_groups) SELECT COUNT(*) FROM ".to_string();

    if met.starts_with("HMDB") {
        query += "metabolites INNER JOIN temp_concat_endogeneity ON metabolites.id = temp_concat_endogeneity.id INNER JOIN temp_concat_functional_groups ON metabolites.id = temp_concat_functional_groups.id INNER JOIN in_tissue ON metabolites.id = in_tissue.id";
    } else {
        query += "lipids";
        return query;
    };
    if matrix == "Positive Mode".to_string() || matrix == "Negative Mode".to_string() {
        query += &format!(
            " WHERE {met_type} {tissue_string}",
            met_type = build_condition_query2(&typ, false).unwrap(),
            tissue_string = tissue_string
        )[..];
        return query;
    } else {
        query += &format!(
            " WHERE {met_type} {functional_group} {tissue_string}",
            met_type = build_condition_query2(&typ, false).unwrap(),
            functional_group = build_condition_query(&adducts, true, "").unwrap(),
            tissue_string = tissue_string
        )[..];
        return query;
    }
}
