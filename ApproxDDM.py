# -*- coding: utf-8 -*-
"""
Implementation of the Approximate Double Description Method

"""
import numpy as np 
from Problem import Problem
from time import process_time as clock




class Node():
    """
    Implements a 3D vector with incidence information         
    """
    node_cnt=0
    def __init__(self):
        self.index=Node.node_cnt
        Node.node_cnt+=1
        self.vector=None
        self.inc=[]
        self.sign=0 # 0: '0'
                    # 1: '-'
                    # 2: '+'
        

    def __str__(self):
        return 'n({0})'.format(self.index)
        
        
    def __repr__(self):
        return 'n({0})'.format(self.index)
        
        
    def info(self):
        print('Node')
        print('  Index          :  {}'.format(self.index))
        print('  Vector         :  {}'.format(self.vector))
        print('  Incidence list :  {}'.format(self.inc))
        print('  Sign           :  {}'.format(V_SIGN[self.sign])) 
        
        
        
        
class Facet():
    """
    Implements "Facet" which corresponds to an inequality and its approximately incident nodes         
    """ 
    def __init__(self):
        self.index=None
        self.node_list=[]
        self.hyperplane=None
        
    def __str__(self):
        return 'f({0})'.format(self.index)
        
        
    def __repr__(self):
        return 'f({0})'.format(self.index)
        
    def info(self):
        print('Facet')
        print('  Index            :  '+self.__str__())
        print('  Node list        :  {}'.format(self.node_list))
        print('  Hyperplane       :  {}'.format(self.hyperplane))  
        
        
        

class ApproxDDM:
    """
    Implements the Approximate Double Description Method         
    """
    def __init__(self, problem):
        self.problem=problem
        self.nodes=[]
        self.facets=[]
        self.matrix=problem.matrix
        self.eps=problem.eps
        self.iter = 0
        self.__init_nodes()
        
        self.time0=0
        self.time1=0
        self.time2=0
        self.time3=0
        self.time4=0
        #self.time5=0
        #self.time6=0
        
        
    def __add_node(self):
        n=Node()
        self.nodes.append(n)
        return n
        
    
    def __remove_node(self,n):
        self.nodes.remove(n)
        del n
        
        
    def __add_facet(self):
        f=Facet()
        self.facets.append(f)
        return f
        
    
    def __remove_facet(self,f):
        self.facets.remove(f)
        del f
        
        
    def __init_nodes(self):
        """
        initialize as simplex
        """
        node_cnt=0
        A=self.matrix[0:4]
        k={0,1,2,3}
        for i in range(4): 
            f=self.__add_facet()
            f.index=i
            f.hyperplane=self.matrix[i]
            n=self.__add_node()         
            n.vector=np.linalg.solve(np.delete(A,i,0),(1+self.eps/2)*np.ones(3))
            n.inc=k.copy()
            n.inc.remove(i)
            self.iter+=1
            
    def __set_v_signs(self):
        hp=self.matrix[self.iter]
        for n in self.nodes:
            #res=hp @ n.vector -> a bit faster:
            res0=hp[0]*n.vector[0]
            res1=hp[1]*n.vector[1]
            res2=hp[2]*n.vector[2]
            res=res0+res1+res2
            if (res < 1):
                n.sign=1 # '-'
            elif (res > 1+self.eps):
                n.sign=2 # '+'
            else:
                n.sign=0 # '0'
                
    def __add_nodes(self):
        hp=self.matrix[self.iter]
        f=self.__add_facet()
        f.index=self.iter
        f.hyperplane=hp
        for w in self.nodes[:]:
            if w.sign==2: # '+'
                for u in self.nodes[:]:
                    if u.sign==1: # '-'
                        if len(set.intersection(u.inc,w.inc))>=2:
                            n=self.__add_node()                        
                            hw=hp @ w.vector 
                            hd=hp @ (u.vector-w.vector)        
                            n.vector = ((1 + self.eps/2 - hw)/ hd) * (u.vector-w.vector) + w.vector
                            n.inc=set.intersection(u.inc,w.inc)
                            n.inc.add(self.iter)
            elif w.sign==0: # '0'
                w.inc.add(self.iter)
                
    
    def __del_nodes(self):
        for w in self.nodes[:]:
            if w.sign==2: # '+'
                self.__remove_node(w)
                
    def step(self):
        if self.iter < self.matrix.shape[0]:
            #print('set_signs')
            t=clock()
            self.__set_v_signs()
            self.time1+=clock()-t
            #print('add_nodes')
            t=clock()
            self.__add_nodes()
            self.time2+=clock()-t
            #print('del_nodes')
            t=clock()
            self.__del_nodes()
            self.time3+=clock()-t
             
            self.iter +=1           
        else:
            print('ApproxVE.step(): no further step to do.')
            
    
    def set_nodelist_of_facets(self):
        for f in self.facets[:]:
            for n in self.nodes:
                if f.index in n.inc:
                    f.node_list.append(n)
            if len(f.node_list)<=2:
                #print('Facet {} removed'.format(f))
                self.__remove_facet(f)
                
                
    def list_facets(self):
        for f in self.facets:
            print('Facet {0} has nodes {1}'.format(f,f.node_list))
            
    
    def list_nodes(self):
        for n in self.nodes:
            print('Node {0} belongs to facet(s) {1}'.format(n,n.facet_list))
    
                
    def update_v_indices(self):
        """
        updates the node indices so that the maximum index
        matches the number of nodes
        """
        ind = 0
        for n in self.nodes:
            n.index=ind
            ind+=1
            
    def run(self):
        t0=clock()
        
        while self.iter < self.matrix.shape[0]:
            #print('ApproxDDM: Processing inequality {}'.format(self.iter))
            self.step()
        self.update_v_indices()
        
        t4=clock()
        self.set_nodelist_of_facets()
        self.time4+=clock()-t4 
        
        self.time0=clock()-t0
        
        
    def kill(self):        
        for n in self.nodes[:]:
            self.__remove_node(n)
        Node.node_cnt=0
        for f in self.facets[:]:
            self.__remove_facet(f)
        Facet.facet_cnt=0
            
    
    def info(self):
        print('--------------------------------------------------')
        print('Method : ApproxDDM')
        print('Problem: {0} with eps= {1}'.format(self.problem.filename,self.eps))
        print('--------------------------------------------------')
        print('CPU-time in seconds .....')
        print('  set_v_signs()  : {}'.format(round(self.time1,3)))
        print('  add_nodes()    : {}'.format(round(self.time2,3)))
        print('  del_nodes()    : {}'.format(round(self.time3,3)))
        print('  get_facets()   : {}'.format(round(self.time4,3)))
        #print('  test           : {}'.format(round(self.time5,3)))
        #print('  test           : {}'.format(round(self.time6,3)))
        print('  --------------')
        print('  total          : {}'.format(round(self.time0,3)))
        print('--------------------------------------------------')
        print('Polytope information .....')
        print('  Nodes   : {}'.format(len(self.nodes)))
        print('  Facets  : {}'.format(len(self.facets)))
        print('--------------------------------------------------')
