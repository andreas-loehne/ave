# -*- coding: utf-8 -*-
"""
Class for Approximate Vertex Enumeration problem instances

"""

import numpy as np 
from scipy.optimize import linprog


class Problem:
    """
    implements an algorithm for Approximate Vertex Enumeration         
    """
    def __init__(self, matrix=None, eps=0.1,filename='myproblem.cvs'):
        self.matrix=matrix
        self.eps=eps
        self.dim=0
        self.filename=filename
        
    
    def __check_matrix(self):
        """
        some trivial tests         
        """
        m=self.matrix.shape[0]
        n=self.matrix.shape[1]
        if m<n+1:
            print('Problem: Matrix invalid: at least dim+1 rows expected')
            exit(1)
        if np.linalg.matrix_rank(self.matrix)!=n:
            print('Problem: Matrix invalid: full rank expected')
            exit(1)
        return m,n

        
    def __get_new_row_alt(self):
        """
        obtain new row by solving the linear program 
        
            max_{t,s} s   s.t.
               A^T t   = 0
                    t >= s e
                e^T t  = 1 
                    s >= 0   
        
            where e = (1,...,1)^T
        """ 
        m=self.matrix.shape[0]
        n=self.matrix.shape[1]
        A=np.zeros((n+1,m+1))
        A[0:n,0:m]=np.transpose(self.matrix)
        A[n,0:m]=1
        B=-np.eye(m+1)
        B[0:m,m]=1
        a=np.zeros((n+1,1))
        a[n]=1
        b=np.zeros((m+1,1))
        c=np.zeros((1,m+1))
        c[0,m]=-1
        res=linprog(c,A_ub=B,b_ub=b,A_eq=A,b_eq=a,method='revised simplex')
        t=res.x[0:m]
        s=res.x[m]
        if s<1e-10:
            print('Polytope is unbounded')
            exit(1)
        # compute new row from lp solution
        new_row = np.zeros(3)
        for i in range(3):  
            new_row = new_row - t[i]*self.matrix[i]
        return new_row
            

    def __get_new_row(self):
        """
        obtain new row by solving the linear program 
        
            max_x c^T   s.t.
               A x <= e
        
            where c is the sum of 3 linearly independent rows

        """ 
        n=self.matrix.shape[1]
        if n==3:
            new_row = self.matrix[0]+self.matrix[1]+self.matrix[2]  
        elif n==2:
            new_row = self.matrix[0]+self.matrix[1]
        else:
            exit(1)             
        m=len(self.matrix)
        mm=min(m,20)       
        res=linprog(c=new_row,A_ub=self.matrix[0:mm],b_ub=np.ones((mm,1)), bounds=(None,None), method='revised simplex')   
        if res.fun == np.NINF:
            print('Polytope is unbounded')
            exit(1)
        new_row=new_row / (res.fun-1) # "minus 1" usually not necessary
        return new_row


    def __prepare_matrix(self):
        """
        matrix is changed: 
            order of rows changed
            one row added
        first 4 rows of resulting matrix describe simplex containing polytope
        exit(1) if polytope is unbounded       
        """
        m, n = self.__check_matrix() 
        if n==3:
            breaker=False
            for i in range(m):
                for j in range(i+1,m):
                    if np.linalg.matrix_rank([self.matrix[i],self.matrix[j],self.matrix[j]])==2:
                        for k in range(j+1,m):
                            if np.linalg.matrix_rank([self.matrix[i],self.matrix[j],self.matrix[k]])==3:
                                breaker=True
                                break
                    if breaker:
                        break
                if breaker:
                    break
            tmp0=self.matrix[i].copy() 
            tmp1=self.matrix[j].copy()
            tmp2=self.matrix[k].copy()
            self.matrix[i]=self.matrix[0]   
            self.matrix[j]=self.matrix[1]
            self.matrix[k]=self.matrix[2]
            self.matrix[0]=tmp0
            self.matrix[1]=tmp1
            self.matrix[2]=tmp2  
            new_row=self.__get_new_row()
            self.matrix=np.append(self.matrix,[new_row],axis=0)
            # change new row with row 4
            tmp=self.matrix[3].copy()
            self.matrix[3]=self.matrix[m].copy()
            self.matrix[m]=tmp.copy()
        elif n==2:
            breaker=False
            for i in range(m):
                for j in range(i+1,m):
                    if np.linalg.matrix_rank([self.matrix[i],self.matrix[j],self.matrix[j]])==2:
                        breaker=True
                        break
                if breaker:
                    break
            tmp0=self.matrix[i].copy() 
            tmp1=self.matrix[j].copy()
            self.matrix[i]=self.matrix[0]   
            self.matrix[j]=self.matrix[1]
            self.matrix[0]=tmp0
            self.matrix[1]=tmp1 
            new_row=self.__get_new_row()
            self.matrix=np.append(self.matrix,[new_row],axis=0)
            # change new row with row 3
            tmp=self.matrix[2].copy()
            self.matrix[2]=self.matrix[m].copy()
            self.matrix[m]=tmp.copy()
            
        else:
            print('matrix must have 2 or 3 columns')
            exit(1)
        
        
    def __load_matrix(self):
        import csv
        with open(self.filename, newline='') as f:
            reader = csv.reader(f,delimiter = ',')    
            A=[]      
            for row in reader:
                dim=len(row)
                A.append([0]*dim)
                for j in range(0,dim):                
                    A[-1][j]=float(row[j])
            self.matrix=np.array(A)
        f.close()
        
    def init(self):
        self.__load_matrix()
        self.__prepare_matrix()
        self.dim=len(self.matrix[0])
        if self.dim > 3 or self.dim < 2:
            print('dimension must be 2 or 3')
            exit(1)
        
        