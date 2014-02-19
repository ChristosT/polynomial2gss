import matrix_coefficients as mc
from sympy.core.symbol import symbols
from sympy.matrices.dense import Matrix, eye, zeros
from sympy.solvers.solvers import solve_linear_system
from sympy import pprint
s=symbols('s')

def is_row_proper(A):
    """
    Test for row properness
    hdr must be full rank! wolovich page 27
    Returns True if A is row proper False otherwise
    
    Example:  TODO
    """
    """
    import matrix_coefficients as mc
    s=symbols('s')
    """
    return mc.highest_row_degree_matrix(A,s).rank() ==min(A.rows,A.cols) # vardu sel 5

def is_column_proper(A):

    return is_row_proper(A.T)


def row_proper(A):
    """
    Reduction of a polynomial matrix to row proper polynomial matrix

    Implementation of the algorithm as shown in page 7 of 
    Vardulakis, A.I.G. (Antonis I.G). Linear Multivariable Control: Algebraic Analysis and Synthesis Methods. 
    Chichester,New York,Brisbane,Toronto,Singapore: John Wiley & Sons, 1991.

    INPUT:
        Polynomial Matrix A
    OUTPUT:
        T:An equivalent matrix in row proper form
        TL the unimodular transformation matrix
    Christos Tsolakis
    
    Example:  TODO
    T=Matrix([[1, s**2, 0], [0, s, 1]])
    is_row_proper(T)
    
    see also Wolovich Linear Multivariable Systems Springer
    """
    T=A.copy()
    Transfomation_Matrix=eye(T.rows)     #initialize transformation matrix
    while  not is_row_proper(T):
        #for i in range(3):
        Thr=mc.highest_row_degree_matrix(T,s)
        #pprint(Thr)
        x=symbols('x0:%d'%(Thr.rows+1))
        Thr=Thr.transpose()
        Thr0=Thr.col_insert(Thr.cols,zeros(Thr.rows,1))
        SOL=solve_linear_system(Thr0,*x)      # solve the system 
        a=Matrix(x)
        D={var:1 for var in x if var not in SOL.keys()}   # put free variables equal to 1
        D.update({var:SOL[var].subs(D) for var in x if var in SOL.keys()})
        a=a.subs(D)
        r0=mc.find_degree(T,s) 
        row_degrees=mc.row_degrees(T,s)                               
        i0=row_degrees.index(max(row_degrees))
        ast=Matrix(1,T.rows,[a[i]*s**(r0-row_degrees[i]) for i in range(T.rows)])
        #pprint (ast)
        TL=eye(T.rows)
        TL[i0*TL.cols]=ast
        #pprint (TL)
        T=TL*T
        T.simplify()
        Transfomation_Matrix=TL*Transfomation_Matrix
        Transfomation_Matrix.simplify()
        #pprint (T)
        #print('----------------------------')
            
    return  Transfomation_Matrix,T

def column_proper(A):
    U,T=row_proper(A.T)
    return U.T,T.T
