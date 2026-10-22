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
    ccs: Option<f64>
}

struct DBNamesIDs {
    origin: String,
    name: String,
    chemicalformula: String,
    id: i64,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct AdductValues {
    pub adjusted_mz: f64,
    pub ccs: Option<f64>,
}
#[derive(Serialize, Debug)]
pub struct ParsedDBData {
    name: String,
    mz: String,
    db_accession: String,
    smiles: String,
    formula: String,
    map: HashMap<String, HashMap<String, AdductValues>>,
    functional_groups: HashMap<String, i64>
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

    // Step 4: Final SQL query
    let sql = format!(
        "SELECT \"{name_col}\" AS matrix_name, ({cases}) AS active_functional_groups FROM user_matrices;",
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

fn functional_group_by_met_id(index: i64, fg_string: &String) -> HashMap<String, i64> {
    let conn: r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager> = get_connection().unwrap();

    let mut stmt = conn.prepare(&format!("SELECT * FROM {fg_string} WHERE rowid = ?1")).unwrap();
    let column_names: Vec<String> = stmt
        .column_names()
        .iter()
        .map(|s| s.to_string())
        .collect();
    let mut rows = stmt.query([index]).unwrap();

    if let Some(row) = rows.next().unwrap() {
        let mut map = HashMap::new();

        for (i, name) in column_names.iter().enumerate() {
            // Try to get the column value as a string
            let value: i64 = row.get(i).unwrap_or_default();

            map.insert(name.to_string(), value);
        }

        map
    } else {
        HashMap::new()
    }
}



pub fn get_db_data(index: usize, origin: String) -> Option<ParsedDBData> {

    let conn: r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager> =
        get_connection().unwrap();

    let matrix_to_groups = get_matrix_functional_groups(&conn).unwrap();
    let adducts = adducts_string(&conn, matrix_to_groups, index, &origin).unwrap();

    let (origin_formula, fg_string, ccs_join) =
        if origin == "metabolites" {
            (
                "metabolites.chemicalformula",
                "functional_groups",
                "LEFT JOIN ccs_predictions AS cp
                    ON cp.metabolite_id = metabolites.id
                   AND cp.adduct = t.adduct",
            )
        } else {
            (
                "lipids.formula",
                "lipids_functional_groups",
                "LEFT JOIN ccs_predictions AS cp ON 1 = 0",
            )
        };

    let metabolite_functional_groups =
        functional_group_by_met_id(index.try_into().unwrap(), &fg_string.to_string());

    let index_i64 = i64::try_from(index).ok()?;

    let query = format!(
        "SELECT
            CAST({origin}.mz AS REAL) +
                CAST(t.deltamass AS REAL) AS adjusted_mz,
            {origin}.name,
            CAST({origin}.mz AS TEXT),
            t.adduct,
            db_accessions.hmdb,
            {origin}.smiles,
            {origin_formula},
            t.mname,
            cp.ccs
         FROM {origin}
         INNER JOIN db_accessions
            ON {origin}.id = db_accessions.id
         INNER JOIN {fg_string}
            ON {origin}.id = {fg_string}.id
         JOIN temp_concat_adducts AS t
            ON t.adduct IN ({adducts})
         {ccs_join}
         WHERE {origin}.id = ?1
         ORDER BY adjusted_mz"
    );

    let mut stmt = conn.prepare_cached(&query).ok()?;

    let rows: Vec<DBData> = stmt
        .query_map([index_i64], |row| {
            Ok(DBData {
                adjusted_mz: row.get(0).unwrap_or_default(),
                name: row.get(1).unwrap_or_default(),
                mz: row.get(2).unwrap_or_default(),
                adduct: row.get(3).unwrap_or_default(),
                db_accession: row.get(4).unwrap_or_default(),
                smiles: row.get(5).unwrap_or_default(),
                chemicalformula: row.get(6).unwrap_or_default(),
                mname: row.get(7).unwrap_or_default(),
                ccs: row.get::<_, Option<f64>>(8).unwrap_or(None),
            })
        })
        .ok()?
        .filter_map(Result::ok)
        .collect();

    rows.first().map(|first| {
        let mut main_hashmap: HashMap<String, HashMap<String, AdductValues>> =
            HashMap::new();

        for row in &rows {
            if row.adjusted_mz == 0.0 {
                continue;
            }

            main_hashmap
                .entry(row.mname.clone())
                .or_default()
                .insert(
                    row.adduct.clone(),
                    AdductValues {
                        adjusted_mz: row.adjusted_mz,
                        ccs: row.ccs,
                    },
                );
        }

        ParsedDBData {
            name: first.name.clone(),
            mz: first.mz.clone(),
            db_accession: first.db_accession.clone(),
            smiles: first.smiles.clone(),
            formula: first.chemicalformula.clone(),
            map: main_hashmap,
            functional_groups: metabolite_functional_groups,
        }
    })
}

pub fn db_ids_and_names(inputvalue: String) -> Vec<(String, String, (String, usize), i64)> {
    let conn = get_connection().unwrap();
    let mut stmt: rusqlite::Statement = conn
        .prepare(
            "SELECT 'metabolites', name, chemicalformula, id FROM metabolites WHERE ((name LIKE ?1) OR (chemicalformula LIKE ?1)) \
                    UNION ALL \
                    SELECT 'lipids', name, formula, id FROM lipids WHERE name LIKE ?1 \
                    LIMIT 20",
        )
        .expect("Query cannot be run");

    let db_iter = stmt
        .query_map([format!("%{}%", inputvalue)], |row: &rusqlite::Row<'_>| {
            Ok(DBNamesIDs {
                origin: row.get(0).unwrap_or("".to_string()),
                name: row.get(1).unwrap_or("".to_string()),
                chemicalformula: row.get(2).unwrap_or("".to_string()),
                id: row.get(3).unwrap_or(0),
            })
        })
        .unwrap();

    let mut map: HashMap<(String, String), (String, usize)> = HashMap::new();

    let matcher = ClangdMatcher::default();

    for (_, item) in db_iter.enumerate() {
        let row: DBNamesIDs = item.unwrap();
        map.insert((row.name, row.chemicalformula), (row.origin, row.id.try_into().unwrap()));
    }

    let mut filtered_names: Vec<(String, String, (String, usize), i64)> = map
        .iter()
        .filter_map(|((name, chemicalformula), id)| {
            let s1 = matcher.fuzzy_match(name, &inputvalue);
            let s2 = matcher.fuzzy_match(chemicalformula, &inputvalue);

            match (s1, s2) {
                (Some(a), Some(b)) => Some((name.clone(), chemicalformula.clone(), id.clone(), a.max(b))),
                (Some(a), None) => Some((name.clone(), chemicalformula.clone(), id.clone(), a)),
                (None, Some(b)) => Some((name.clone(), chemicalformula.clone(), id.clone(), b)),
                (None, None) => None,
            }
        })
        .collect();

    // Sort first by fuzzy match score (descending), then by length, then alphabetically
    filtered_names.sort_by(|a, b| {
        b.3.cmp(&a.3)                   // fuzzy score (largest first)
            .then_with(|| a.0.len().cmp(&b.0.len())) // shorter names first
            .then_with(|| a.0.cmp(&b.0))             // alphabetical
    });

    filtered_names
}
