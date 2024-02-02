import tkinter as tk

def main():
    def on_element_click(element_name):
        current_text = search_bar.get()
        new_text = current_text + " " + element_name
        search_bar.delete(0, tk.END)
        search_bar.insert(0, new_text)

    root = tk.Tk()
    root.title("DB GUI demo")
    root.geometry("800x600")

    top_frame = tk.Frame(root)
    top_frame.pack(side='top')
    #top_frame.grid(row=0, column=0, columnspan=18)  # Adjust as needed
    search_bar = tk.Entry(top_frame)
    search_bar.pack(pady=10)

    bottom_frame = tk.Frame(root)
    bottom_frame.pack(fill='x',expand=True)

    for element in elements:
        btn = tk.Button(bottom_frame, text=element["symbol"], 
                        command=lambda e=element["name"]: on_element_click(e))
        btn.grid(row=element["row"], column=element["column"], sticky="nsew")

    for i in range(1, 19):  # Configure columns
        bottom_frame.grid_columnconfigure(i, weight=1)

    for i in range(1, 8):  # Configure rows
        bottom_frame.grid_rowconfigure(i, weight=1)

    root.mainloop()

    return 0

#

elements = [
    {"symbol": "H", "name": "Hydrogen", "row": 1, "column": 1},
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
    {"symbol": "Xe", "name": "Xenon", "row": 5, "column": 18}
    # ... Continue for other elements
]

#
if __name__ == '__main__':
  main()
