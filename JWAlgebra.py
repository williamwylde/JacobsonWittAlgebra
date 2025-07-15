from sage.all import *
from itertools import *

#Formats element of the JW Algebra as a variable name
def Wn_format(alph, N):
    return "x{}D{}x".format(alph,N).replace('(','H').replace(')', 'b').replace(',', 'c').replace(' ', '')

#Translates variable name into readable string
def Wn_read(elem):
    elem_string = str(elem)
    return elem_string.replace('H', '(').replace('b',')').replace('c',',').replace('x', '')

#Returns a dictionary of an element of the JW Algebra
def Wn_basis_element(alph, N):
    return {'name': Wn_format(alph, N), 'alpha':alph, 'i':N}

#Returns a dictionary of Wn basis elements given a prime 'p' and rank 'n'
def Wn_basis(p, n):
    A = list(product(range(p), repeat=n)) #gives n-tuples of 0<=a_i<p
    basis = {}
    for N in range(n):
        for alph in A:
            basis_element = Wn_basis_element(alph, N)
            basis[basis_element['name']]=basis_element
    return basis

#Function to add alphas in the basis element, and check there is no overflow modulo p
def Wn_basis_alpha_add_helper(alphas, p):
    output = tuple(map(sum, zip(*alphas)))
    if any(k < 0 or k >= p for k in output):
        return None
    else:
        return output

#Construct [(alpha)Di, (beta)Dj] according to a basis and prime p
def Wn_bracket(a, b, basis, p):
    n = len(a['alpha'])
    i = a['i']
    j = b['i']
    alpha = a['alpha']
    beta = b['alpha']
    w_i = tuple([0 if k!=i else -1 for k in range(0,n)])
    w_j = tuple([0 if k!=j else -1 for k in range(0,n)])
    ab_wj = Wn_basis_alpha_add_helper([alpha,beta,w_j], p)
    ab_wi = Wn_basis_alpha_add_helper([alpha,beta,w_i], p)
    result = {}
    if ab_wj is not None:
        br_1 = Wn_basis_element(ab_wj, i)
        result[br_1['name']] = -1*alpha[j]
    if ab_wi is not None:
        br_2 = Wn_basis_element(ab_wi, j)
        if br_2['name'] in result:
            result[br_2['name']] += beta[i]
        else:
            result[br_2['name']] = beta[i]
    return result

#Generates bracket dictionary readable by SAGE
def Wn_generate_bracket(W_basis, p):
    bracket_output = {}
    for A in W_basis:
        for B in W_basis:
            if str(A) < str(B):
                br_output = Wn_bracket(W_basis[A],W_basis[B],W_basis,p)
                if br_output != {}:
                    bracket_output = bracket_output | {(W_basis[A]['name'],W_basis[B]['name']):br_output}
    return bracket_output

#Construct a Jacobson Witt algebra W_n over a field Z_p as a SAGE LieAlgebra object
def JacobsonWitt(p, n):
    Wn = Wn_basis(p, n)
    bracket = Wn_generate_bracket(Wn,p)
    #print(bracket)
    K = GF(p)
    return LieAlgebra(K, list(Wn.keys()), bracket)

def JW(*args):
    return JacobsonWitt(*args)

#ad function - very useful
def ad(x, y):
    return x*y - y*x

#Writes a file of generators of the center of the UEA of the (p, n) JW Algebra
def Write_UEA_central_generators(filename, p, n, ngens):
    L = JacobsonWitt(p, n)
    U = L.pbw_basis()
    Z = U.center()
    g = Z.algebra_generators()
    with open(filename, "a") as f:
        f.write('\nGenerators of Z(U(W{})) in Z{}\n'.format(n, p))
        for i in range(0,ngens):
            f.write('{}th generator is: '.format(str(i)) + Wn_read(U(g[i])) + '\n')

#Absolutely horrible code to turn a PBW monomial into a dictionary for use with LaTeX code further down
def Strip_PBW_monomial(mon):
    if "^" in mon:
        exp = int(mon.split("^")[1]) ##Add error checking
    else:
        exp = None
    if "*" in mon:
        coeff = int(mon.split("*")[0])
    else:
        coeff = None
    #get alpha#
    temp1 = mon.split("PBW['")[1]
    temp2 = temp1.split("D")[0]
    alpha = eval(temp2)
    #get indices#
    temp1 = mon.split("D")[1]
    temp2 = temp1.split("']")[0]
    i = int(temp2)
    return {'coeff': coeff, 'i': i, 'exp': exp, 'alpha': alpha}

#Returns the LaTeX code for a given PBW monomial
def LaTeX_PBW_mon(mon):
    coeff = mon['coeff']
    ind = mon['i']
    alpha = mon['alpha']
    exp = mon['exp']
    if coeff == None:
        coeff = ""
    if exp == None:
        exp = ""
    return " $[{co}\\delta_{{{ind}}}^{{{al}}}]^{{{ex}}}$ ".format(co = coeff, ind=ind, al=alpha, ex=exp)

#Splits an element of the UEA by '*' and '+'
def Split_PBW_by_operator(PBW):
    output = []
    plus = PBW.split(" ")
    for word in plus:
        if word.count("PBW") <=1:
            output.append(word)
        else:
            subwords = word.split("*P")
            for j in range(1, len(subwords)):
                subwords[j] = "P" + subwords[j]
            split_subword = ['*'] * (len(subwords) * 2 - 1)
            split_subword[0::2] = subwords
            for thing in split_subword:
                output.append(thing)
    return output

#Returns the LaTeX code for an element of the UEA
def LaTeX_centraliser(elem):
    file = Split_PBW_by_operator(elem)
    output = ""
    for word in file:
        if "PBW" in word:
            output += LaTeX_PBW_mon(Strip_PBW_monomial(word))
        else:
            output += word
    return output

#Writes LaTeX file with generators of Z(U(W_n)), halt stops the program if the next generator is algebraically dependent on the previous generators
def Write_generators_as_LaTeX(filename, p, n, ngens, halt=True):
    L = JacobsonWitt(p, n)
    U = L.pbw_basis()
    Z = U.center()
    gen = Z.algebra_generators()
    with open(filename, "a") as f:
        f.write('\\section*{{Generators of Z(U(W{})) in Z{}}}\n'.format(n, p))
        for i in range(0,ngens):
            g = gen[i]
            g_read = Wn_read(U(g))
            if ("*" in str(g) or "^" in str(g)) and halt: #g looks like Z[0], Z[1], etc. Once we see Z[0]^2 or Z[0]*Z[1] we have exhausted the algebraically independent generators of Z(U)
                break
            f.write('{}th generator is: {}, written as:'.format(str(i), str(g).replace("^", "**")) + g_read.replace("*", " * ").replace("^", "**") + '\\\\')
            f.write('$$' + LaTeX_centraliser(g_read).replace("$", "") + '$$\\\\ \n')