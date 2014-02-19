from sympy import *
import reduction 
import matrix_coefficients as mc
from itertools import product
s=symbols('s')
def ALGO21(As,Bs,Cs,Ds,do_test):
    
    #-------------------STEP 1------------------------------
    r=As.rows
    p=Cs.rows
    m=Ds.cols

    Ts=BlockMatrix([[As,Bs,zeros(As.rows,p)],
                    [-Cs,Ds,eye(p)],
                    [zeros(m,Cs.cols),-eye(m),zeros(m,p)]]).as_mutable()
                  
    rt=r+p+m
    T_over_I=BlockMatrix([[Ts],[eye(rt)]]).as_mutable()
    H=reduction.column_proper(T_over_I)
    Qs_over_Rs=BlockMatrix([[H[1][0:rt,0:rt]],[H[0]]]).as_mutable()    
    #-------------------END STEP 1------------------------------
    
    #-------------------STEP 2------------------------------
    qi=mc.column_degrees(Qs_over_Rs,s)
    Ss=BlockDiagMatrix(*[Matrix(q+1,1,[s**(q-k) for k in range(q+1)]) for q in qi]).as_mutable() ## POS??POSS?   
    
    #find Qc
    numvar=rt*Ss.rows
    a=symbols('a0:%d'%numvar)
    Qc=Matrix(rt,Ss.rows,a)
    RHS_matrix=(Qc*Ss).as_mutable()
    Qs=Qs_over_Rs[0:rt,0:rt]

    SOL=dict()
    for i,j in  product(range(Qs.rows),range(Qs.cols)):                 
        SOL.update(solve_undetermined_coeffs(Eq(Qs[i,j],RHS_matrix[i,j]),a,s))
    Qc=Qc.subs(SOL) 

    #find Rc
    numvar=rt*Ss.rows
    Rc=Matrix(rt,Ss.rows,a)
    RHS_matrix=(Rc*Ss).as_mutable()
    Rs=H[0]
    SOL=dict()
    for i,j in  product(range(Rs.rows),range(Rs.cols)):                 
        SOL.update(solve_undetermined_coeffs(Eq(Rs[i,j],RHS_matrix[i,j]),a,s))
    Rc=Rc.subs(SOL) 
    #-------------------END STEP 2------------------------------
    
    #-------------------STEP 3------------------------------

    def E0block(i,j):
        if i==j and i!=0:
            return 1
        else:
            return 0
    
    E0=BlockDiagMatrix(*[Matrix(q+1,q+1, lambda i,j: E0block(i,j)) for q in qi]).as_mutable()
    A0=BlockDiagMatrix(*[Matrix(q+1,q+1, lambda i,j: KroneckerDelta(i-1,j)) for q in qi]).as_mutable()
    B0=BlockDiagMatrix(*[Matrix(q+1,1, lambda i,j: KroneckerDelta(i,j)) for q in qi]).as_mutable()
    n=sum(qi)+rt
    C0=eye(n)
    V=BlockMatrix([[zeros(p,r),zeros(p,m),eye(p)]]).as_mutable()
    U=BlockMatrix([[zeros(r,m)],[zeros(p,m)],[eye(m)]]).as_mutable()
    
    E=SparseMatrix(E0)
    A=SparseMatrix(A0-B0*Qc)
    B=SparseMatrix(B0*U)
    C=SparseMatrix(V*Rc)
    #E=E0.as_mutable()
    #pprint(B0*Qc)
    #A=(simplify(A0-B0*Qc)).as_mutable()
    #B=(B0*U).as_mutable()
    #C=(V*Rc).as_mutable()
     #-------------------END STEP 3------------------------------

    if do_test==True:
        test_result=test(B0,Ts,E,A,B,C,U,V,Rc,Rs,Ss,p,m,rt)
    else:
        test_result= 'not done' 
    
    D=zeros(p,m) 
    return E,A,B,C,D,test_result

def test(B0,Ts,E,A,B,C,U,V,Rc,Rs,Ss,p,m,rt):
       
    D1=BlockDiagMatrix(B0,eye(p)).as_mutable()
    D2=BlockMatrix([[Ts,U],[-V,zeros(p,m)]]).as_mutable()
    D3=BlockMatrix([[s*E-A,B],[-C,zeros(C.rows,B.cols)]]).as_mutable()
    D4Block11=simplify((Ss*Rs.inv())).as_mutable()
    D4=BlockDiagMatrix(D4Block11,eye(m)).as_mutable()
    return simplify(D1*D2)==simplify(D3*D4) 
    
    






