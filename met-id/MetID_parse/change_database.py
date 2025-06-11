import csv
import sqlite3



file_path = 'names_identifiers.tsv'
data = {}

with open(file_path, 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    for row in reader:
        print(row)
        try:
            key, value = row[0], row[1]
            data[key] = value
        except IndexError:
            print(row[0])


conn = sqlite3.connect('msms_db.db')
cursor = conn.cursor()

# Add the new column 'foz'
#cursor.execute('ALTER TABLE MSMS ADD COLUMN identifier TEXT')


# Select the values from the 'foo' column
cursor.execute('SELECT id, name FROM MSMS')
rows = cursor.fetchall()

# Update the 'foz' column based on the dictionary translation
for row in rows:
    rowid = row[0]
    foo_value = row[1]
    foz_value = data.get(foo_value, None) # Translate the value using the dictionary
    if foz_value is not None:
        cursor.execute('UPDATE MSMS SET identifier = ? WHERE id = ?', (foz_value, rowid))

# Commit the changes
conn.commit()

# Close the connection
conn.close()