import tkinter as tk
from tkinter import ttk
import sqlite3

def main():
    app = DBApp()
    app.run()
    
#
class DBApp:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("DB GUI demo")
        self.root.geometry("800x600")
        self.setup_search_UI()
        
    def run(self):
        self.root.mainloop()
    
    def setup_search_UI(self):
        # Search bar
        search_frame = tk.Frame(self.root)
        search_frame.pack(side='top')
        #top_frame.grid(row=0, column=0, columnspan=18)  # Adjust as needed
        self.search_var = tk.StringVar()
        self.search_bar = tk.Entry(search_frame, textvariable=self.search_var)
        self.search_bar.pack(side='left', pady=10, expand=True)
        self.search_button = tk.Button(search_frame, text="Search", command=self.perform_search)
        #self.search_bar.bind("<Return>", self.perform_search)
        self.search_button.pack(side='left', padx=10, pady=10)

        # Periodic table
        periodic_frame = tk.Frame(self.root)
        periodic_frame.pack(side='bottom', fill='x',expand=True)
        elements = [{"symbol": "H", "name": "Hydrogen", "row": 1, "column": 1},
                    {"symbol": "He", "name": "Helium", "row": 1, "column": 18},
                    {"symbol": "Li", "name": "Lithium", "row": 2, "column": 1},
                    {"symbol": "Be", "name": "Beryllium", "row": 2, "column": 2},
                    {"symbol": "B", "name": "Boron", "row": 2, "column": 13},
                    {"symbol": "C", "name": "Carbon", "row": 2, "column": 14},
                    {"symbol": "N", "name": "Nitrogen", "row": 2, "column": 15},
                    {"symbol": "O", "name": "Oxygen", "row": 2, "column": 16},
                    {"symbol": "F", "name": "Fluorine", "row": 2, "column": 17},
                    {"symbol": "Ne", "name": "Neon", "row": 2, "column": 18},
                    {"symbol": "Na", "name": "Sodium", "row": 3, "column": 1},
                    {"symbol": "Mg", "name": "Magnesium", "row": 3, "column": 2},
                    {"symbol": "Al", "name": "Aluminum", "row": 3, "column": 13},
                    {"symbol": "Si", "name": "Silicon", "row": 3, "column": 14},
                    {"symbol": "P", "name": "Phosphorus", "row": 3, "column": 15},
                    {"symbol": "S", "name": "Sulfur", "row": 3, "column": 16},
                    {"symbol": "Cl", "name": "Chlorine", "row": 3, "column": 17},
                    {"symbol": "Ar", "name": "Argon", "row": 3, "column": 18},
                    {"symbol": "K", "name": "Potassium", "row": 4, "column": 1},
                    {"symbol": "Ca", "name": "Calcium", "row": 4, "column": 2},
                    {"symbol": "Sc", "name": "Scandium", "row": 4, "column": 3},
                    {"symbol": "Ti", "name": "Titanium", "row": 4, "column": 4},
                    {"symbol": "V", "name": "Vanadium", "row": 4, "column": 5},
                    {"symbol": "Cr", "name": "Chromium", "row": 4, "column": 6},
                    {"symbol": "Mn", "name": "Manganese", "row": 4, "column": 7},
                    {"symbol": "Fe", "name": "Iron", "row": 4, "column": 8},
                    {"symbol": "Co", "name": "Cobalt", "row": 4, "column": 9},
                    {"symbol": "Ni", "name": "Nickel", "row": 4, "column": 10},
                    {"symbol": "Cu", "name": "Copper", "row": 4, "column": 11},
                    {"symbol": "Zn", "name": "Zinc", "row": 4, "column": 12},
                    {"symbol": "Ga", "name": "Gallium", "row": 4, "column": 13},
                    {"symbol": "Ge", "name": "Germanium", "row": 4, "column": 14},
                    {"symbol": "As", "name": "Arsenic", "row": 4, "column": 15},
                    {"symbol": "Se", "name": "Selenium", "row": 4, "column": 16},
                    {"symbol": "Br", "name": "Bromine", "row": 4, "column": 17},
                    {"symbol": "Kr", "name": "Krypton", "row": 4, "column": 18},
                    {"symbol": "Rb", "name": "Rubidium", "row": 5, "column": 1},
                    {"symbol": "Sr", "name": "Strontium", "row": 5, "column": 2},
                    {"symbol": "Y", "name": "Yttrium", "row": 5, "column": 3},
                    {"symbol": "Zr", "name": "Zirconium", "row": 5, "column": 4},
                    {"symbol": "Nb", "name": "Niobium", "row": 5, "column": 5},
                    {"symbol": "Mo", "name": "Molybdenum", "row": 5, "column": 6},
                    {"symbol": "Tc", "name": "Technetium", "row": 5, "column": 7},
                    {"symbol": "Ru", "name": "Ruthenium", "row": 5, "column": 8},
                    {"symbol": "Rh", "name": "Rhodium", "row": 5, "column": 9},
                    {"symbol": "Pd", "name": "Palladium", "row": 5, "column": 10},
                    {"symbol": "Ag", "name": "Silver", "row": 5, "column": 11},
                    {"symbol": "Cd", "name": "Cadmium", "row": 5, "column": 12},
                    {"symbol": "In", "name": "Indium", "row": 5, "column": 13},
                    {"symbol": "Sn", "name": "Tin", "row": 5, "column": 14},
                    {"symbol": "Sb", "name": "Antimony", "row": 5, "column": 15},
                    {"symbol": "Te", "name": "Tellurium", "row": 5, "column": 16},
                    {"symbol": "I", "name": "Iodine", "row": 5, "column": 17},
                    {"symbol": "Xe", "name": "Xenon", "row": 5, "column": 18},
                    {"symbol": "Cs", "name": "Cesium", "row": 6, "column": 1},
                    {"symbol": "Ba", "name": "Barium", "row": 6, "column": 2},
                    # Lanthanides (57-71) are skipped
                    {"symbol": "Hf", "name": "Hafnium", "row": 6, "column": 4},
                    {"symbol": "Ta", "name": "Tantalum", "row": 6, "column": 5},
                    {"symbol": "W", "name": "Tungsten", "row": 6, "column": 6},
                    {"symbol": "Re", "name": "Rhenium", "row": 6, "column": 7},
                    {"symbol": "Os", "name": "Osmium", "row": 6, "column": 8},
                    {"symbol": "Ir", "name": "Iridium", "row": 6, "column": 9},
                    {"symbol": "Pt", "name": "Platinum", "row": 6, "column": 10},
                    {"symbol": "Au", "name": "Gold", "row": 6, "column": 11},
                    {"symbol": "Hg", "name": "Mercury", "row": 6, "column": 12},
                    {"symbol": "Tl", "name": "Thallium", "row": 6, "column": 13},
                    {"symbol": "Pb", "name": "Lead", "row": 6, "column": 14},
                    {"symbol": "Bi", "name": "Bismuth", "row": 6, "column": 15},
                    {"symbol": "Po", "name": "Polonium", "row": 6, "column": 16},
                    {"symbol": "At", "name": "Astatine", "row": 6, "column": 17},
                    {"symbol": "Rn", "name": "Radon", "row": 6, "column": 18},
                    {"symbol": "Fr", "name": "Francium", "row": 7, "column": 1},
                    {"symbol": "Ra", "name": "Radium", "row": 7, "column": 2},
                    # Actinides (89-103) are skipped
                    {"symbol": "Rf", "name": "Rutherfordium", "row": 7, "column": 4},
                    {"symbol": "Db", "name": "Dubnium", "row": 7, "column": 5},
                    {"symbol": "Sg", "name": "Seaborgium", "row": 7, "column": 6},
                    {"symbol": "Bh", "name": "Bohrium", "row": 7, "column": 7},
                    {"symbol": "Hs", "name": "Hassium", "row": 7, "column": 8},
                    {"symbol": "Mt", "name": "Meitnerium", "row": 7, "column": 9},
                    {"symbol": "Ds", "name": "Darmstadtium", "row": 7, "column": 10},
                    {"symbol": "Rg", "name": "Roentgenium", "row": 7, "column": 11},
                    {"symbol": "Cn", "name": "Copernicium", "row": 7, "column": 12},
                    {"symbol": "Nh", "name": "Nihonium", "row": 7, "column": 13},
                    {"symbol": "Fl", "name": "Flerovium", "row": 7, "column": 14},
                    {"symbol": "Mc", "name": "Moscovium", "row": 7, "column": 15},
                    {"symbol": "Lv", "name": "Livermorium", "row": 7, "column": 16},
                    {"symbol": "Ts", "name": "Tennessine", "row": 7, "column": 17},
                    {"symbol": "Og", "name": "Oganesson", "row": 7, "column": 18}]
        elements += [{"symbol": "La", "name": "Lanthanum", "row": 9, "column": 3},
                     {"symbol": "Ce", "name": "Cerium", "row": 9, "column": 4},
                     {"symbol": "Pr", "name": "Praseodymium", "row": 9, "column": 5},
                     {"symbol": "Nd", "name": "Neodymium", "row": 9, "column": 6},
                     {"symbol": "Pm", "name": "Promethium", "row": 9, "column": 7},
                     {"symbol": "Sm", "name": "Samarium", "row": 9, "column": 8},
                     {"symbol": "Eu", "name": "Europium", "row": 9, "column": 9},
                     {"symbol": "Gd", "name": "Gadolinium", "row": 9, "column": 10},
                     {"symbol": "Tb", "name": "Terbium", "row": 9, "column": 11},
                     {"symbol": "Dy", "name": "Dysprosium", "row": 9, "column": 12},
                     {"symbol": "Ho", "name": "Holmium", "row": 9, "column": 13},
                     {"symbol": "Er", "name": "Erbium", "row": 9, "column": 14},
                     {"symbol": "Tm", "name": "Thulium", "row": 9, "column": 15},
                     {"symbol": "Yb", "name": "Ytterbium", "row": 9, "column": 16},
                     {"symbol": "Lu", "name": "Lutetium", "row": 9, "column": 17}]
        elements += [{"symbol": "Ac", "name": "Actinium", "row": 10, "column": 3},
                     {"symbol": "Th", "name": "Thorium", "row": 10, "column": 4},
                     {"symbol": "Pa", "name": "Protactinium", "row": 10, "column": 5},
                     {"symbol": "U", "name": "Uranium", "row": 10, "column": 6},
                     {"symbol": "Np", "name": "Neptunium", "row": 10, "column": 7},
                     {"symbol": "Pu", "name": "Plutonium", "row": 10, "column": 8},
                     {"symbol": "Am", "name": "Americium", "row": 10, "column": 9},
                     {"symbol": "Cm", "name": "Curium", "row": 10, "column": 10},
                     {"symbol": "Bk", "name": "Berkelium", "row": 10, "column": 11},
                     {"symbol": "Cf", "name": "Californium", "row": 10, "column": 12},
                     {"symbol": "Es", "name": "Einsteinium", "row": 10, "column": 13},
                     {"symbol": "Fm", "name": "Fermium", "row": 10, "column": 14},
                     {"symbol": "Md", "name": "Mendelevium", "row": 10, "column": 15},
                     {"symbol": "No", "name": "Nobelium", "row": 10, "column": 16},
                     {"symbol": "Lr", "name": "Lawrencium", "row": 10, "column": 17}]

        def on_element_click(input_str):
            current_text = self.search_bar.get()
            new_text = current_text + " " + input_str
            self.search_bar.delete(0, tk.END)
            self.search_bar.insert(0, new_text)   
        for element in elements:
            btn = tk.Button(periodic_frame, text=element["symbol"], 
                            command=lambda e=element["symbol"]: on_element_click(e))
            btn.grid(row=element["row"], column=element["column"], sticky="news")
        for i in range(1, 19):  # Configure columns
            periodic_frame.grid_columnconfigure(i, weight=1)
        for i in range(1, 11):  # Configure rows
            periodic_frame.grid_rowconfigure(i, weight=1)
    
    def perform_search(self):#, event=None):
        search_term = self.search_var.get()
        self.show_search_results(search_term)
    
    def show_search_results(self, search_term):
        # Create a popup window
        popup = tk.Toplevel(self.root)
        popup.title("Search Results")
        popup.geometry("800x600")
        
        # Treeview widget
        tree = ttk.Treeview(popup, columns=("Project Name", "Type", "File Path", "File Length"), show="headings")
        tree.heading("Project Name", text="Project Name", command=lambda: self.treeview_sort_column(tree, "Project Name", False))
        tree.heading("Type", text="Type", command=lambda: self.treeview_sort_column(tree, "Type", False))
        tree.heading("File Path", text="File Path", command=lambda: self.treeview_sort_column(tree, "File Path", False))
        tree.heading("File Length", text="File Length", command=lambda: self.treeview_sort_column(tree, "File Length", False))
        tree.column("Project Name", anchor=tk.W, width=100)
        tree.column("Type", anchor=tk.W, width=100)
        tree.column("File Path", anchor=tk.W, width=200)
        tree.column("File Length", anchor=tk.W, width=100)
        tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)        
        scrollbar = ttk.Scrollbar(popup, orient=tk.VERTICAL, command=tree.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        tree.configure(yscrollcommand=scrollbar.set)

        # Query the database and insert data into the Treeview
        #load_data()
        conn = sqlite3.connect('my.db')
        cur = conn.cursor()
        cur.execute("SELECT project_name, type, file_path, file_length FROM project_files")
        rows = cur.fetchall()
        for row in rows:
            tree.insert("", tk.END, values=row)
        export_button = tk.Button(popup, text="Export", command=self.export_selected(tree,cur))
        export_button.pack()
        conn.close()

    def treeview_sort_column(self, treeview, col, reverse):
        l = [(treeview.set(k, col), k) for k in treeview.get_children('')]
        l.sort(reverse=reverse)
        for index, (val, k) in enumerate(l):
            treeview.move(k, '', index)
        treeview.heading(col, command=lambda _col=col: self.treeview_sort_column(treeview, _col, not reverse))

    def export_selected(self, treeview, cur):
        selected_items = treeview.selection()
        data_to_export = []
        for item_id in selected_items:
            item = treeview.item(item_id)
            values = item['values']
            cur.execute("SELECT * FROM your_table WHERE id=?", (values[0],))  # Adjust query as needed
            data_to_export.append(cur.fetchone())
        # Write the data to a text file
        with open('exported_data.txt', 'w') as f:
            for data in data_to_export:
                f.write(f"{data}\n")
        print("Data exported successfully to 'exported_data.txt'.")

#
if __name__ == '__main__':
  main()