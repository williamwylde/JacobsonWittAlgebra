import tkinter as tk
import pyperclip
from JWAlgebra import *

class ControlInterface():
    def __init__(self, master):
        self.master = master
        self.master.title("Jacobson Witt Algebra")
        self.construct_JW_generator_interface()
        self.construct_basis_element_interface()

    def update(self):
        self.master.update_idletasks()
        self.master.update()
    
    def construct_JW_generator_interface(self):
        tk.Label(self.master, text="n: ").grid(row=0,column=0)
        self.set_n = tk.Entry(self.master, width=2)
        self.set_n.grid(row=0, column=1)

        tk.Label(self.master, text="p: ").grid(row=1,column=0)
        self.set_p = tk.Entry(self.master, width=2)
        self.set_p.grid(row=1, column=1)

        tk.Button(self.master, text='Construct', command=self.construct_JW_generator).grid(row=2, column=1)
        tk.Label(self.master, text="Generate W_n over Z_p").grid(row=3,column=1)
        tk.Label(self.master, text="L = W_n, U = U(L), Z = Z(U(L))").grid(row=4,column=1)

    def construct_JW_generator(self):
        cmd = "L = JW({},{}); U = L.pbw_basis(); Z = U.center(); Zgen = Z.algebra_generators()".format(self.set_p.get(), self.set_n.get())
        print("Executing: " + cmd)
        exec(cmd, globals())
        self.update()

    def construct_basis_element_interface(self):
        tk.Label(self.master, text="alpha: ").grid(row=0,column=2)
        self.set_alpha = tk.Entry(self.master, width=10)
        self.set_alpha.grid(row=0, column=3)

        tk.Label(self.master, text="i: ").grid(row=1,column=2)
        self.set_i = tk.Entry(self.master, width=2)
        self.set_i.grid(row=1, column=3)
        tk.Label(self.master, text="Copy generator in: ").grid(row=2,column=2)
        tk.Button(self.master, text='L', command=self.construct_basis_element_in_L).grid(row=2, column=3)
        tk.Button(self.master, text='U', command=self.construct_basis_element_in_U).grid(row=3, column=3)
    def construct_basis_element_in_L(self):
        a = eval(self.set_alpha.get())
        i = int(self.set_i.get())
        cmd = "L.basis()['{}']".format(Wn_format(a, i))
        pyperclip.copy(cmd)
        self.update()
    
    def construct_basis_element_in_U(self):
        a = eval(self.set_alpha.get())
        i = int(self.set_i.get())
        cmd = "U(L.basis()['{}'])".format(Wn_format(a, i))
        pyperclip.copy(cmd)
        self.update()


#Leave this at the bottom
#Set global variables
global L, U, Z, Zgen
L = JW(2,2); U = L.pbw_basis(); Z = U.center(); Zgen = Z.algebra_generators()

root = tk.Tk()
CI = ControlInterface(root)