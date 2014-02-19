from sympy import *
import reduction 
import matrix_coefficients as mc
from itertools import product
s=symbols('s')
def ALGO24(As,Bs,Cs,Ds,do_test):
    
    #-------------------STEP 1------------------------------
    r=As.rows
    p=Cs.rows
    m=Ds.cols

    Ts=BlockMatrix([[As,Bs,zeros(As.rows,p)],
                    [-Cs,Ds,eye(p)],
                    [zeros(m,Cs.cols),-eye(m),zeros(m,p)]]).as_mutable()
    
    rt=r+p+m
    I_over_T=BlockMatrix([[eye(rt)],[Ts]]).as_mutable()
    H=reduction.column_proper(I_over_T)
    Qs_over_Rs=H[1] 

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
    Rs=Qs_over_Rs[rt:2*rt,0:rt]
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
    
    E=BlockDiagMatrix(E0,zeros(rt,rt)).as_mutable()
    A=BlockMatrix([[A0-B0*Qc,-B0],[Rc,zeros(Rc.rows,B0.cols)]]).as_mutable()
    B=BlockMatrix([[zeros(n,m)],[U]]).as_mutable()
    C=BlockMatrix([[zeros(p,n),V]]).as_mutable()

    #-------------------END STEP 3------------------------------
    
    if do_test==True:
        test_result=test(B0,Ts,E,E0,A0,A,B,C,U,V,Rc,Rs,Qs,Ss,p,m,rt,n)
    else:
        test_result= 'not done' 
    
    D=zeros(p,m) 
    return E,A,B,C,D,test_result
            
    
def test(B0,Ts,E,E0,A0,A,B,C,U,V,Rc,Rs,Qs,Ss,p,m,rt,n):

    D1=BlockMatrix([[zeros(n,rt+p)],[eye(rt +p)]]).as_mutable()
    D2=BlockMatrix([[Ts,U],[-V,zeros(p,m)]]).as_mutable()
    D3=BlockMatrix([[s*E-A,B],[-C,zeros(C.rows,B.cols)]]).as_mutable()
    D4Block1=simplify((-Ss*(Qs.inv()))).as_mutable()
    D4Block1=D4Block1.col_join(eye(rt))
    D4=BlockDiagMatrix(D4Block1,eye(m)).as_mutable()
    

    return simplify(D1*D2) == simplify(D3*D4)
    

