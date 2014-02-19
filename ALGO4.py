from sympy.core.symbol import symbols
from sympy.matrices import Matrix, diag, SparseMatrix, eye
from sympy.matrices.expressions.blockmatrix import BlockDiagMatrix ,BlockMatrix
from sympy.core.function  import expand
from sympy.functions.special.tensor_functions import KroneckerDelta 
from sympy.core.relational import Eq
from sympy.polys.polytools import poly,div
from sympy.solvers.solvers import solve , solve_undetermined_coeffs
from sympy.simplify.simplify import simplify
from sympy import numer,denom
from itertools import product

from reduction import *
import matrix_coefficients as mc
s=symbols('s')
def ALGO4(As,Bs,Cs,Ds,do_test):
    
    
    
#-------------------STEP 1------------------------------
    if not is_row_proper(As):
        Us,Ast=row_proper(As)
        
    else:
        Us,Ast=eye(As.rows),As
    
    Bst=Us*Bs
    Bst=expand(Bst)
    r=Ast.cols
    #-------------------STEP 2------------------------------
    K=simplify(Ast.inv()*Bst)  #very important 
    Ys=zeros(K.shape)
    for i,j in  product(range(K.rows),range(K.cols)):
        Ys[i,j],q=div(numer(K[i,j]),denom(K[i,j]))       
    
    B_hat=Bst-Ast*Ys
    
    #-------------------END STEP 2------------------------------
    #-------------------STEP 3------------------------------
    Psi=diag(*[[s**( mc.row_degrees(Ast,s)[j]  -i -1) for i in range( mc.row_degrees(Ast,s)[j])] for j in range(r)]).T
    S=diag(*[s**(rho) for rho in mc.row_degrees(Ast,s)])
    Ahr=mc.highest_row_degree_matrix(Ast,s)
    Help=Ast-S*Ahr
    
    SOL={}
    numvar=Psi.rows*Psi.cols
    alr=symbols('a0:%d'%numvar)
    Alr=Matrix(Psi.cols,Psi.rows,alr)
    RHS=Psi*Alr
    for i,j in  product(range(Help.rows),range(Help.cols)):                 #diagonal explain later
        SOL.update(solve_undetermined_coeffs(Eq(Help[i,j],RHS[i,j]),alr,s))
    
    Alr=Alr.subs(SOL)    #substitute(SOL)
    
    
    Aoc=Matrix(BlockDiagMatrix(*[Matrix(rho, rho, lambda i,j: KroneckerDelta(i+1,j))for rho in mc.row_degrees(Ast,s)]))
    Boc=eye(sum(mc.row_degrees(Ast,s)))
    Coc=Matrix(BlockDiagMatrix(*[SparseMatrix(1,rho,{(x,0):1 if x==0 else 0 for x in range(rho)}) for rho in mc.row_degrees(Ast,s)]))

    A0=Aoc-Alr*Ahr.inv()*Matrix(Coc)
    C0=Ahr.inv()*Coc



    SOL={}
    numvar=Psi.cols*Bst.cols
    b0=symbols('b0:%d'%numvar)
    B0=Matrix(Psi.cols,Bst.cols,b0)
    RHS=Psi*B0

    for i,j in  product(range(B_hat.rows),range(B_hat.cols)):                 #diagonal explain later
        SOL.update(solve_undetermined_coeffs(Eq(B_hat[i,j],RHS[i,j]),b0,s))
    B0=B0.subs(SOL)    #substitute(SOL)

    LHS_matrix=simplify(Cs*C0)                                        #left hand side of the equation (1)

    sI_A=s*eye(A0.cols)- A0
    max_degree=mc.find_degree(LHS_matrix,s)                   #get the degree of the matrix at the LHS
                                                  #which is also the maximum degree for the coefficients of Λ(s)
    #---------------------------Creating Matrices Λ(s) and C -------------------------------------

    Lamda=[]
    numvar=((max_degree))*A0.cols
    a=symbols('a0:%d'%numvar)
    for i in range(A0.cols):                                        # paratirisi den douleuei to prin giat;i otra oxi diagonios
        p=sum(a[n +i*(max_degree)]*s**n for n in range(max_degree))  # we want variables one degree lower because we are multiplying by first order monomials
        Lamda.append(p)
        
    Lamda=Matrix(Cs.rows,A0.cols,Lamda)                                #convert the list to Matrix
    
    c=symbols('c0:%d'%(Lamda.rows*Lamda.cols))
    C=Matrix(Lamda.rows,Lamda.cols,c)
    #-----------------------------------------
    
    RHS_matrix=Lamda*sI_A +C                            #right hand side of the equation (1)
     
    '''
    -----------Converting equation (1) to a system of linear -----------
    -----------equations, comparing the coefficients of the  -----------
    -----------polynomials in both sides of the equation (1) -----------
    '''
     EQ=[Eq(LHS_matrix[i,j],expand(RHS_matrix[i,j])) for i,j in product(range(LHS_matrix.rows),range(LHS_matrix.cols)) ]
    eq=[]
    
    for equation in EQ:
        RHS=poly((equation).rhs,s).all_coeffs() #simplify necessary ?
        LHS=poly((equation).lhs,s).all_coeffs()
        if len(RHS)>len(LHS):                      # we add zero for each missing coefficient (greater than  the degree of LHS)
            LHS=(len(RHS)-len(LHS))*[0] + LHS
        eq=eq+[Eq(LHS[i],RHS[i]) for i in range(len(RHS))]
    
    SOL=solve(eq,a+c)   # the coefficients of Λ and C
    
    
    #----------substitute the solution in the matrices----------------
    C=C.subs(SOL)     #substitute(SOL)
    Lamda=Lamda.subs(SOL)    #substitute(SOL)
    Js=(Ds+Cs*Ys+Lamda*(B0));
    # for output compatibility
    E=eye(A0.cols)
    if do_test==True:
        test_result=test(As,Bs,Cs,Ds,Us,sI_A,Js,Ys,B0,C0,C,Psi,Lamda)
    else:
        test_result= 'not done' 
    return (E,A0,B0,C,Js,test_result)


def test(As,Bs,Cs,Ds,Us,sI_A,Js,Ys,B0,C0,C,Psi,Lamda):

    D1=BlockMatrix([[Us.inv()*Psi,zeros(Us.rows,1)],[-Lamda,eye(Lamda.rows)]]).as_mutable()
    D2=BlockMatrix([[sI_A,B0],[-C,Js]]).as_mutable()
    D3=BlockMatrix([[As,Bs],[-Cs,Ds]]).as_mutable()
    D4=BlockMatrix([[C0,-Ys],[zeros(Ys.cols,C0.cols),eye(Ys.cols)]]).as_mutable()
    return expand(simplify(D1*D2))==expand(D3*D4)



"""
s=symbols('s')
As=Matrix(2,2,[s+1 ,s,s**2,s**4 ]);
#TODO include tests
r=As.rows
Bs=Matrix([[s+1],[s**5+2*s**2+s+3]]); # to 5 htan 3
Cs=Matrix([[-s**2-3*s-1,1]]);
Ds=Matrix([[s**3 +s +2]]);

s=symbols('s')
As=Matrix(2,2,[s+1 ,s**3+2*s**2,s**2 +3*s+2,s**4+4*s**3 +4*s**2 +s+2]);
#TODO include tests
r=As.rows
Bs=Matrix([[s**2+1],[s**3+2*s**2+s+3]]);
Cs=Matrix([[-s**2-3*s-1,-s**4 -4*s**3-4*s**2+1]]);
Ds=Matrix([[s**3+2*s**2 +s +2]]);

s=symbols('s')
As=Matrix(3,3,[s+1 ,s,s**2,s**4, 5,s+6,s**2 +s**3,9*s**2,10 ]);
#TODO include tests
r=As.rows
Bs=Matrix([[s+1],[s**5+2*s**2+s+3],[s*2 +s**3 +1]]); # to 5 htan 3
Cs=Matrix([[-s**2-3*s-1,1,s**2 -s]]);
Ds=Matrix([[s**3 +s +2]]);
"""
