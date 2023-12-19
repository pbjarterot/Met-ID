import sqlite3

# Replace 'your_database.db' with the path to your SQLite database file
database_path = 'msms_db.db'

# SQL command to alter the table
sql_command = """
ALTER TABLE MSMS
ADD COLUMN matrix TEXT DEFAULT 'FMP-10';
"""

# Connect to the SQLite database
conn = sqlite3.connect(database_path)

try:
    # Create a cursor object and execute the SQL command
    cursor = conn.cursor()
    cursor.execute(sql_command)

    # Commit the changes
    conn.commit()
    print("Column added successfully.")

except sqlite3.Error as error:
    print("Error while executing SQL script:", error)

finally:
    # Close the connection
    if conn:
        conn.close()