from sympy import Matrix

def check_input(A,B,C,D):
    r=A.rows
    m=B.cols
    p=C.rows
    if A.is_square and A.rank()==r and B.rows==r and C.cols==r and D.rows==p and D.cols==m:
        return True
    else:
        return False
    
