import sqlite3

conn = sqlite3.connect('db.db')
cursor = conn.cursor()

# Query all employees
cursor.execute("""
SELECT (CAST(metabolites.mz AS REAL) + 
CASE WHEN adducts.adduct IN ('M+FMP10', 'M+2FMP10a', 'M+2FMP10b', 'M+3FMP10a', 'M+3FMP10b', 'M+3FMP10c') 

THEN CAST(adducts.deltamass AS REAL) ELSE 0 END) AS adjusted_mz,metabolites.name, adducts.adduct, db_accessions.hmdb, metabolites.smiles, metabolites.chemicalformula 
FROM metabolites 
INNER JOIN db_accessions ON metabolites.id = db_accessions.id 
INNER JOIN endogeneity ON metabolites.id = endogeneity.id 
INNER JOIN functional_groups ON metabolites.id = functional_groups.id 
INNER JOIN derivatized_by ON metabolites.id = derivatized_by.id 
CROSS JOIN (SELECT 'M+FMP10' AS mname 
UNION SELECT 'M+2FMP10a' AS mname
UNION SELECT 'M+2FMP10b' 
UNION SELECT 'M+3FMP10a' 
UNION SELECT 'M+3FMP10b' 
UNION SELECT 'M+3FMP10c') 
AS m LEFT JOIN adducts ON m.mname = adducts.adduct 

WHERE CASE WHEN adducts.adduct='M+FMP10' 
  AND derivatized_by.fmp >= 1
  AND (endogeneity.endogenous = 1 OR endogeneity.exogenous = 1 OR endogeneity.unspecified = 1) 
  AND (functional_groups.phenols + functional_groups.primaryamines >= 1 ) 
  THEN 1 
WHEN adducts.adduct IN ('M+FMP10', 'M+2FMP10a', 'M+2FMP10b') 
  AND derivatized_by.fmp >= 2
  AND (endogeneity.endogenous = 1 OR endogeneity.exogenous = 1 OR endogeneity.unspecified = 1) 
  AND (functional_groups.phenols + functional_groups.primaryamines >= 2) 
  THEN 1 

WHEN adducts.adduct IN ('M+FMP10', 'M+2FMP10a', 'M+2FMP10b', 'M+3FMP10a', 'M+3FMP10b', 'M+3FMP10c') 
  AND derivatized_by.fmp > 2 
  AND (endogeneity.endogenous = 1 OR endogeneity.exogenous = 1 OR endogeneity.unspecified = 1) 
  AND (functional_groups.phenols + functional_groups.primaryamines > 2) 
  THEN 1 
  ELSE 0 
  END=1 
ORDER BY CAST(metabolites.mz AS REAL) + CAST(adducts.deltamass AS REAL)




""")
rows = cursor.fetchall()

# Print the results
for row in rows:
    if row[0] < 413.3 and row[0] > 413.2:
      print(row)

# Close the connection
conn.close()











"""
SELECT (CAST(metabolites.mz AS REAL) + 
CASE WHEN adducts.adduct IN ('M+FMP10', 'M+2FMP10a', 'M+2FMP10b', 'M+3FMP10a', 'M+3FMP10b', 'M+3FMP10c') 

THEN CAST(adducts.deltamass AS REAL) ELSE 0 END) AS adjusted_mz,metabolites.name, adducts.adduct, db_accessions.hmdb, metabolites.smiles, metabolites.chemicalformula 
FROM metabolites 
INNER JOIN db_accessions ON metabolites.id = db_accessions.id 
INNER JOIN endogeneity ON metabolites.id = endogeneity.id 
INNER JOIN functional_groups ON metabolites.id = functional_groups.id 
INNER JOIN derivatized_by ON metabolites.id = derivatized_by.id 
CROSS JOIN (SELECT 'M+FMP10' AS mname 
UNION SELECT 'M+2FMP10a' AS mname
UNION SELECT 'M+2FMP10b' 
UNION SELECT 'M+3FMP10a' 
UNION SELECT 'M+3FMP10b' 
UNION SELECT 'M+3FMP10c') 
AS m LEFT JOIN adducts ON m.mname = adducts.adduct 

WHERE CASE WHEN adducts.adduct='M+FMP10' 
  AND derivatized_by.fmp = 1 
  AND (endogeneity.endogenous = 1 OR endogeneity.exogenous = 1 OR endogeneity.unspecified = 1) 
  AND (functional_groups.phenols = 1 OR functional_groups.primaryamines = 1) 
  THEN 1 
WHEN adducts.adduct IN ('M+FMP10', 'M+2FMP10a', M+2FMP10b') 
  AND derivatized_by.fmp = 2 
  AND CAST(metabolites.mz AS REAL) + CAST(adducts.deltamass AS REAL) < 1200 
  AND CAST(metabolites.mz AS REAL) + CAST(adducts.deltamass AS REAL) > 0  
  AND (endogeneity.endogenous = 1 OR endogeneity.exogenous = 1 OR endogeneity.unspecified = 1) 
  AND (functional_groups.phenols = 1 OR functional_groups.primaryamines = 1) 
  THEN 1 

WHEN adducts.adduct IN ('M+FMP10', 'M+2FMP10a', 'M+2FMP10b', 'M+3FMP10a', 'M+3FMP10b', 'M+3FMP10c') 
  AND derivatized_by.fmp>2 
  AND CAST(metabolites.mz AS REAL) + CAST(adducts.deltamass AS REAL) < 1200 
  AND CAST(metabolites.mz AS REAL) + CAST(adducts.deltamass AS REAL) > 0 
  AND (endogeneity.endogenous = 1 OR endogeneity.exogenous = 1 OR endogeneity.unspecified = 1) 
  AND (functional_groups.phenols = 1 OR functional_groups.primaryamines = 1) 
  THEN 1 
  ELSE 0 
  END=1 
ORDER BY CAST(metabolites.mz AS REAL) + CAST(adducts.deltamass AS REAL)




"""