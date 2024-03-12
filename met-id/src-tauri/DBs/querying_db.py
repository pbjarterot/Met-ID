import sqlite3

# Replace 'your_database.db' with the path to your SQLite database file
database_path = 'msms_db.db'

# SQL command to alter the table
sql_command = """
SELECT name, adduct FROM MSMS;
"""

# Connect to the SQLite database
conn = sqlite3.connect(database_path)

try:
    # Create a cursor object and execute the SQL command
    cursor = conn.cursor()
    cursor.execute(sql_command)

    names = []
    for i in cursor.fetchall():
        if i[1] in ["1A", "2A", "2B"]:
            names.append(f"{i[0]}, {i[1]}")

    print(len(list(set(names))))


except sqlite3.Error as error:
    print("Error while executing SQL script:", error)

finally:
    # Close the connection
    if conn:
        conn.close()