import tkinter as tk
from tkinter import ttk
import sqlite3

class PeriodicTableApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Periodic Table Search")
        self.root.geometry("800x600")
        
        self.setup_ui()
    
    def setup_ui(self):
        # Search bar
        self.search_var = tk.StringVar()
        self.search_bar = tk.Entry(self.root, textvariable=self.search_var)
        self.search_bar.pack(pady=10)
        
        # Search button
        self.search_button = tk.Button(self.root, text="Search", command=self.perform_search)
        self.search_button.pack(pady=5)
    
    def perform_search(self):
        search_term = self.search_var.get()
        self.show_search_results(search_term)
    
    def show_search_results(self, search_term):
        # Create a popup window
        popup = tk.Toplevel(self.root)
        popup.title("Search Results")
        popup.geometry("600x400")
        
        # Treeview widget
        tree = ttk.Treeview(popup, columns=("Project Name", "Type", "File Path", "File Length"), show="headings")
        tree.heading("Project Name", text="Project Name")
        tree.heading("Type", text="Type")
        tree.heading("File Path", text="File Path")
        tree.heading("File Length", text="File Length")
        
        tree.pack(expand=True, fill="both")
        
        # Query the database and insert data into the Treeview
        conn = sqlite3.connect('my.db')
        cur = conn.cursor()
        query = "SELECT project_name, type, file_path, file_length FROM project_files WHERE project_name LIKE ?"
        cur.execute(query, ('%' + search_term + '%',))
        
        for row in cur.fetchall():
            tree.insert("", tk.END, values=row)
        
        conn.close()

root = tk.Tk()
app = PeriodicTableApp(root)
root.mainloop()
