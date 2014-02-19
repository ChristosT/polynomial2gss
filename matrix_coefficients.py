from sympy.core.symbol import symbols
from sympy.matrices import Matrix, SparseMatrix , BlockMatrix
from sympy.polys.polytools import Poly, degree
from sympy.core.function  import expand


def find_degree(M,s):
    '''
    Computes the degree of a polynomial matrix
    M : polynomial Matrix
    s: variable of the polynomial matrix
    
    Example:  TODO
    s=symbols('s')
    T=Matrix([[1, s**2, 0], [0, s, 1]])
    find_degree(T)
    '''
    return max(degree(poly,s) for poly in M if poly!=0 ) # for 0?
    
def row_degrees(A,s):
    """
    return a list of the degrees of the row of matrix
    """
    return [max(degree(poly,s)  for poly in A[i,:] if poly!=0 ) for i in range(A.rows)]

def column_degrees(A,s):
    return row_degrees(A.T,s)

def full_coeffs(p,D,s):
    '''
    Computes coefficients of a polynomial matrix
    p : polynomial
    D :degree  up to which we want to complete the list 	with zeroes
    s: variable of the polynomial matrix
    Example:  TODO
    '''
    p=Poly(p,s,domain='QQ')  # in order to use the all_coeffs method of polynomial class
    c=p.all_coeffs()
    if len(c)==D+1:
        return c
    else:
        difference=D+1-len(c)
        return difference*[0] +c


def matrix_coeffs(A,s):
    '''
    Returns the coefficients of a polynomial matrix
    including those which the coefficient is the Zero matrix
    A : polynomial Matrix
    D :degree of matrix up to which we want to complete the list with zeroes
    Example:  TODO
    '''
    return [Matrix(A.rows,A.cols,[full_coeffs(a,find_degree(A,s),s)[i] for a in A]) for i in range(0,find_degree(A,s)+1)]

def highest_row_degree_matrix(T,s):
    """
    Returns a matrix containing with the coefficients of the a polynomial 
    including those which the coefficient is the Zero matrix
    
    computes the highest row degree coefficient matrix of poynomial matrix A
    A : polynomial Matrix
    D :degree of matrix up to which we want to complete the list with zeroes
    Example:  TODO
    """
    s=symbols('s')
    RD=row_degrees(T,s)
    return  Matrix(T.rows,T.cols, lambda i,j: expand(T[i,j]).coeff(s,RD[i]))

def lowest_row_degree_matrix(A,s):
    """
    Returns a matrix containing with the coefficients of the a polynomial 
    including those which the coefficient is the Zero matrix
    
    computes the lowest row degree coefficient matrix of poynomial matrix A
    A : polynomial Matrix
    D :degree of matrix up to which we want to complete the list with zeroes
    Example:  TODO
    """
    H=Matrix(A.rows,A.cols,[Poly(p,s) for p in A])
    RD=row_degrees(H,s)
    #lowest_row_degrees=[min(min(poly.as_dict().keys())  for poly in H[i,:] if poly!=0 ) for i in range (H.rows)]
    
    return  Matrix(BlockMatrix(H.rows,H.cols, lambda i,j: Matrix(full_coeffs(H[i,j],RD[i],s))))


