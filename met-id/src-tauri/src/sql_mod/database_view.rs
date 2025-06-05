use crate::database_setup::get_connection;
use fuzzy_matcher::{clangd::ClangdMatcher, FuzzyMatcher};
use rusqlite::Result;
use serde::Serialize;
use std::collections::HashMap;

#[derive(Debug)]
struct DBData {
    adjusted_mz: f64,
    mz: String,
    name: String,
    adduct: String,
    db_accession: String,
    smiles: String,
    chemicalformula: String,
    mname: String,
}

struct DBNamesIDs {
    origin: String,
    name: String,
    id: usize,
}
#[derive(Serialize, Debug)]
pub struct ParsedDBData {
    name: String,
    mz: String,
    db_accession: String,
    smiles: String,
    formula: String,
    map: HashMap<String, HashMap<String, f64>>,
}

fn get_matrix_functional_groups(
    conn: &r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager>,
) -> Result<HashMap<String, String>> {
    // Step 1: Get column names
    let mut stmt = conn.prepare("PRAGMA table_info(matrices);")?;
    let columns: Vec<String> = stmt
        .query_map([], |row| Ok(row.get::<_, String>(1)?))?
        .filter_map(|res| res.ok())
        .collect();

    // Step 2: Separate the matrix name column (first) and functional group columns
    let (matrix_name_col, functional_group_cols) = columns.split_first().expect("No columns found");

    // Step 3: Build CASE expressions for each FG column
    let case_expr = functional_group_cols
        .iter()
        .map(|col| {
            format!("CASE WHEN IFNULL(CAST(\"{col}\" AS INTEGER), 0) = 1 THEN '{col},' ELSE '' END")
        })
        .collect::<Vec<_>>()
        .join(" || ");

    // Step 4: Final SQL query
    let sql = format!(
        "SELECT \"{name_col}\" AS matrix_name, ({cases}) AS active_functional_groups FROM matrices;",
        name_col = matrix_name_col,
        cases = case_expr
    );

    let mut stmt = conn.prepare(&sql)?;

    // Step 5: Collect results into a HashMap
    let rows = stmt.query_map([], |row| {
        let matrix_name: String = row.get("matrix_name")?;
        let raw_groups: String = row.get("active_functional_groups")?;
        Ok((matrix_name, raw_groups.trim_end_matches(',').to_string()))
    })?;

    let mut map = HashMap::new();
    for row in rows {
        let (name, groups) = row?;
        map.insert(name, groups);
    }
    map.remove("Positive Mode");
    map.remove("Negative Mode");

    Ok(map)
}

fn adducts_string(
    conn: &r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager>,
    map: HashMap<String, String>,
    id: usize,
    origin: &String,
) -> Result<String> {
    let mut conn2: r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager> =
        get_connection().unwrap();
    let mut select_stmt = "".to_string();

    for (idx, (matrix, fgs)) in map.iter().enumerate() {
        select_stmt += &format!("(mname = \"{matrix}\" AND ");
        let mut some_string = "SELECT ".to_string();
        let split: Vec<&str> = fgs.split(",").collect();

        for (idx, part) in split.iter().enumerate() {
            if idx == split.len() - 1 {
                some_string += &format!("\"{part}\" ");
            } else {
                some_string += &format!("\"{part}\" + ");
            }
        }
        if origin == &"metabolites".to_string() {
            some_string += &format!("FROM functional_groups WHERE id = {id}");
        } else {
            some_string += &format!("FROM lipids_functional_groups WHERE id = {id}");
        }

        select_stmt += &format!("({some_string}) >= temp_concat_adducts.numfunctionalgroups)");

        if !(idx == map.len() - 1) {
            select_stmt += " OR "
        }
    }

    let tx: rusqlite::Transaction<'_> = conn2.transaction().unwrap();
    tx.execute("DROP TABLE IF EXISTS temp_concat_adducts", [])
        .unwrap();
    tx.execute("CREATE TEMP TABLE temp_concat_adducts AS SELECT * FROM adducts UNION ALL SELECT * FROM user_adducts", []).unwrap();

    let sql = format!(
        "
        SELECT adduct
        FROM temp_concat_adducts
        WHERE (temp_concat_adducts.mname IN (\"Positive Mode\", \"Negative Mode\"))
        OR
        {select_stmt}
        "
    );

    let mut stmt = conn.prepare(&sql)?;

    let rows = stmt.query_map([], |row| {
        let adduct: String = row.get(0)?;
        Ok(adduct)
    })?;
    let mut adducts_string = String::new();
    for row in rows {
        adducts_string += &format!("\"{}\", ", row.unwrap());
    }

    //remove the last space and comma
    adducts_string.pop();
    adducts_string.pop();

    Ok(adducts_string)
}

pub fn get_db_data(index: usize, origin: String) -> Option<ParsedDBData> {
    let conn: r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager> =
        get_connection().unwrap();

    //get all unique matrix names and their adducts and numfunctionalgroups
    let matrix_to_groups: HashMap<String, String> = get_matrix_functional_groups(&conn).unwrap();
    let adducts: String = adducts_string(&conn, matrix_to_groups, index, &origin).unwrap();
    let origin_formula: String;

    if origin == "metabolites".to_string() {
        origin_formula = "metabolites.chemicalformula".to_string();
    } else {
        origin_formula = "lipids.formula".to_string();
    }

    let fg_string: String;
    if origin == "metabolites".to_string() {
        fg_string = "functional_groups".to_string();
    } else {
        fg_string = "lipids_functional_groups".to_string();
    }

    let query: String = format!("SELECT (CAST({origin}.mz AS REAL) + CASE WHEN temp_concat_adducts.adduct IN ({adducts}) THEN CAST(temp_concat_adducts.deltamass AS REAL) END) AS adjusted_mz, \
                                {origin}.name, CAST({origin}.mz AS TEXT), temp_concat_adducts.adduct, db_accessions.hmdb, {origin}.smiles, {origin_formula}, temp_concat_adducts.mname FROM {origin} \
                                INNER JOIN db_accessions ON {origin}.id = db_accessions.id \
                                INNER JOIN {fg_string} ON {origin}.id = {fg_string}.id \
                                JOIN temp_concat_adducts ON 1=1 \
                                WHERE {origin}.id = {index} ORDER BY adjusted_mz");

    let mut stmt: rusqlite::Statement = conn.prepare(&query).expect("Query cannot be run");
    let rows: Vec<DBData> = stmt
        .query_map([], |row| {
            Ok(DBData {
                adjusted_mz: row.get(0).unwrap_or_default(),
                name: row.get(1).unwrap_or_default(),
                mz: row.get(2).unwrap_or_default(),
                adduct: row.get(3).unwrap_or_default(),
                db_accession: row.get(4).unwrap_or_default(),
                smiles: row.get(5).unwrap_or_default(),
                chemicalformula: row.get(6).unwrap_or_default(),
                mname: row.get(7).unwrap_or_default(),
            })
        })
        .unwrap()
        .filter_map(Result::ok)
        .collect();

    rows.first().map(|first| {
        let mut main_hashmap: HashMap<String, HashMap<String, f64>> = HashMap::new();

        for row in rows.iter().skip(1) {
            if row.adjusted_mz == 0.0 {
                continue;
            }

            main_hashmap
                .entry(row.mname.clone())
                .or_default()
                .insert(row.adduct.clone(), row.adjusted_mz);
        }
        ParsedDBData {
            name: first.name.clone(),
            mz: first.mz.clone(),
            db_accession: first.db_accession.clone(),
            smiles: first.smiles.clone(),
            formula: first.chemicalformula.clone(),
            map: main_hashmap,
        }
    })
}

pub fn db_ids_and_names(inputvalue: String) -> Vec<(String, (String, usize), i64)> {
    let conn = get_connection().unwrap();
    let mut stmt: rusqlite::Statement = conn
        .prepare(
            "SELECT 'metabolites', name, id FROM metabolites WHERE name LIKE ?1  \
                        UNION ALL \
                        SELECT 'lipids', name, id FROM lipids WHERE name LIKE ?1 LIMIT 20",
        )
        .expect("Query cannot be run");
    let db_iter = stmt
        .query_map([format!("%{}%", inputvalue)], |row: &rusqlite::Row<'_>| {
            Ok(DBNamesIDs {
                origin: row.get(0).unwrap_or("".to_string()),
                name: row.get(1).unwrap_or("".to_string()),
                id: row.get(2).unwrap_or(0),
            })
        })
        .unwrap();
    let mut map: HashMap<String, (String, usize)> = HashMap::new();

    let matcher = ClangdMatcher::default();

    for (_, item) in db_iter.enumerate() {
        let row: DBNamesIDs = item.unwrap();
        map.insert(row.name, (row.origin, row.id));
    }

    // Filter and rank names based on fuzzy match score
    let mut filtered_names: Vec<(String, (String, usize), i64)> = map
        .iter()
        .filter_map(|(name, id)| {
            matcher
                .fuzzy_match(name, &inputvalue)
                .map(|score| (name.clone(), id.clone(), score))
        })
        .collect();

    // Sort first by fuzzy match score (descending), then by length, then alphabetically
    filtered_names.sort_by(|a, b| {
        b.2.cmp(&a.2)
            .then_with(|| a.0.len().cmp(&b.0.len()))
            .then_with(|| a.0.cmp(&b.0))
    });

    filtered_names
}
