import sqlite3
import os

# Path to your raw file
raw_file_path = './db_list.txt'

# Connect to SQLite database (or create it if it doesn't exist)
conn = sqlite3.connect('my.db')
cur = conn.cursor()

# Create a table
cur.execute('''
CREATE TABLE IF NOT EXISTS project_files (
    project_name TEXT,
    type TEXT,
    file_path TEXT,
    file_length INTEGER
)
''')

# Function to calculate file length (number of lines here, can be modified to file size)
def file_length(file_path):
    try:
        with open(file_path, 'r') as f:
            return len(f.readlines())
    except FileNotFoundError:
        return 0

# Parse the raw file and insert data into the database
current_project = ''
current_type = ''
with open(raw_file_path, 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith('#'):  # Project name line
            current_project = ' '.join(line.split()[1:])
        elif line.startswith('['):  # Type line
            current_type = line[1:-1].strip()
        elif line:  # File path line
            path = line
            length = file_length(path)  # Calculate file length
            # Insert into the database
            cur.execute('INSERT INTO project_files (project_name, type, file_path, file_length) VALUES (?, ?, ?, ?)',
                        (current_project, current_type, path, length))

# Commit changes and close the connection
conn.commit()
conn.close()

print("Database created and populated successfully.")
