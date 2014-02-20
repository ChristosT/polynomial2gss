from sympy import *
import numpy as np

import matrix_coefficients as mc
s=symbols('s')
def ALGO11(As,Bs,Cs,Ds,do_test):
    #=======================STEP 1======================================
    r=As.rows
    m=Bs.cols
    p=Cs.rows
    
    Ps=BlockMatrix([[As,Bs],[-Cs,Ds]]).as_mutable()
    List_P=(mc.matrix_coeffs(Ps,s))
    List_P.reverse()
    q=len(List_P)
    k=List_P[0].rows
    l=List_P[0].cols
    Zero_Mat=zeros(List_P[0].rows,List_P[0].cols)  #list of zeroes
    #==================END==STEP 1======================================
    
    #=======================STEP 2======================================
    #Haskel matrices #wikipedia ???
    
    #PiE
    A=[]
    for i in range(2,q):
        A.append([List_P[j] for j in range(i,q)])
    
    for i in range(1,len(A)):
        A[i]=A[i] +[Zero_Mat]*i
        
    PiE=BlockMatrix(A).as_mutable()

    A=[]
    for i in range(3,q+1):
        A.append([List_P[j] for j in range(i,q)])
        A[i-3]=A[i-3] +[Zero_Mat]*(i-2)
    
    PiA=BlockMatrix(A)
    PiB=BlockMatrix((q-2),1,List_P[2:])
    PiC=BlockMatrix(1,q-2,List_P[2:])
    
    rE=PiE.rank()
    #===================END=====STEP 2==================================
    #=======================STEP 3======================================
    J=PiE.rref()[1]                 #pivot columns
    I=(PiE.transpose()).rref()[1]   #pivot rows
    
    PE=Matrix(np.mat(PiE)[np.ix_(I,J)])
    PA=Matrix(np.mat(PiA)[np.ix_(I,J)])
    PB=Matrix(np.mat(PiB)[np.ix_(I,range(l))])
    PC=Matrix(np.mat(PiC)[np.ix_(range(k),J)])
    #===================END=STEP 3======================================
    #=======================STEP 4======================================
    Lambda=r+rE +p +m
    E=BlockMatrix([[List_P[1],PC],[PB,PA],[zeros(Lambda-k-PB.rows,PB.cols),zeros(Lambda-k-PB.rows,PA.cols)]]).as_mutable()
    E=E.row_join(zeros(Lambda,Lambda-l-PC.cols))
    A=BlockMatrix([[-List_P[0],zeros(k,PE.cols),zeros(k-p,Lambda-l-PE.cols).col_join(-eye(p))],
    [zeros(PE.rows,l),PE,zeros(PE.rows,p)],
    [zeros(m,l-m).row_join(eye(m)),zeros(m,PE.cols),zeros(m,p)]]).as_mutable()
    
    
    B=BlockMatrix([[zeros(r+p+rE,m)],[eye(m)]]).as_mutable()
    C=BlockMatrix([[zeros(p,r+m+rE),eye(p)]]).as_mutable()
    D=zeros(p,m)
    #===============END==STEP 4=========================================

    if do_test==True:
        test_result=test(rE,r,p,m,PE,PA,PC,As,Bs,Cs,Ds,E,A,B,C,D)
    else:
        test_result= 'not done' 
    return E,A,B,C,D,test_result
    
def test(rE,r,p,m,PE,PA,PC,As,Bs,Cs,Ds,E,A,B,C,D):
    D1Block2=simplify(PC*s*(PE-s*PA).inv())
    D1Block3=BlockMatrix([[Bs],[Ds]]).as_mutable()
    D1Block4=BlockMatrix([[zeros(r,p)],[eye(p)]]).as_mutable()
    D1=BlockMatrix([[eye(r+p),D1Block2,D1Block3,D1Block4]]).as_mutable()
    D2=BlockMatrix(2,2,[s*E-A,B,-C,D]).as_mutable()

    D3=BlockMatrix(2,2,[As,Bs,-Cs,Ds]).as_mutable()
    D4=BlockMatrix([[eye(r),zeros(r,rE+2*p),zeros(r,m)],
                    [zeros(m,r),zeros(m,rE+2*p),eye(m)]]).as_mutable()
    return simplify(D1*D2)==simplify(D3*D4)





