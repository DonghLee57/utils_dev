import tkinter as tk
from tkinter import ttk
import sqlite3

class DatabaseApp:
    def __init__(self, root):
        self.root = root
        self.root.title("SQLite3 Database Viewer")
        self.root.geometry("600x400")
        self.create_widgets()

    def create_widgets(self):
        self.tree = ttk.Treeview(self.root, columns=("Project Name", "Type", "File Path", "File Length"), show="headings")
        self.tree.heading("Project Name", text="Project Name", command=lambda: self.treeview_sort_column(self.tree, "Project Name", False))
        self.tree.heading("Type", text="Type", command=lambda: self.treeview_sort_column(self.tree, "Type", False))
        self.tree.heading("File Path", text="File Path", command=lambda: self.treeview_sort_column(self.tree, "File Path", False))
        self.tree.heading("File Length", text="File Length", command=lambda: self.treeview_sort_column(self.tree, "File Length", False))

        self.tree.column("Project Name", anchor=tk.W, width=100)
        self.tree.column("Type", anchor=tk.W, width=100)
        self.tree.column("File Path", anchor=tk.W, width=200)
        self.tree.column("File Length", anchor=tk.W, width=100)

        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.scrollbar = ttk.Scrollbar(self.root, orient=tk.VERTICAL, command=self.tree.yview)
        self.scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.tree.configure(yscrollcommand=self.scrollbar.set)

        self.load_data()

    def load_data(self):
        conn = sqlite3.connect('my.db')
        cur = conn.cursor()
        cur.execute("SELECT project_name, type, file_path, file_length FROM project_files")
        rows = cur.fetchall()

        for row in rows:
            self.tree.insert("", tk.END, values=row)

        conn.close()

    def treeview_sort_column(self, tv, col, reverse):
        l = [(tv.set(k, col), k) for k in tv.get_children('')]
        l.sort(reverse=reverse)

        for index, (val, k) in enumerate(l):
            tv.move(k, '', index)

        tv.heading(col, command=lambda _col=col: self.treeview_sort_column(tv, _col, not reverse))

root = tk.Tk()
app = DatabaseApp(root)
root.mainloop()
