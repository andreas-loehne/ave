# -*- coding: utf-8 -*-
"""
Implementation of a Graph-based Approximate VE Algorithm

"""
import numpy as np 
from sys import exit
from time import process_time as clock
from time import sleep
from random import random
from random import seed
from inspect import currentframe, getframeinfo
from Graph import HalfEdge
from Graph import Vertex
from Graph import Face
from Graph import Graph
from Graph import H_TYPE
from Graph import V_SIGN
from Problem import Problem


class GraphAlg: 
    
    def __init__(self, problem):
        self.problem=problem
        self.graph  = None
        self.matrix=problem.matrix
        self.dim=problem.dim
        self.eps=problem.eps
        self.pertubation1=0
        self.pertubation2=0
        self.__init_graph()
        self.iter = 0
        
        self.time0=0
        self.time1=0
        self.time2=0
        self.time3=0
        self.time4=0
        self.time5=0
        self.time6=0
        self.time7=0
        self.time8=0
        self.time9=0
        

    def __init_graph(self):
        """
        initialize as simplex
        """
        self.graph=Graph(self.dim)
        A=self.matrix[0:self.dim+1]
        for i in range(self.dim+1):           
             self.graph.vertices[i].vector=np.linalg.solve(np.delete(A,i,0),(1+self.eps/2)*np.ones(self.dim))
        self.iter=self.dim+1
        
        
    def __set_v_signs(self):
        hp=self.matrix[self.iter]
        for v in self.graph.vertices:
            #res=hp @ v.vector #-> a bit faster:
            if self.dim==3:
                res0=hp[0]*v.vector[0]
                res1=hp[1]*v.vector[1]
                res2=hp[2]*v.vector[2]
                res=res0+res1+res2
            else:
                res0=hp[0]*v.vector[0]
                res1=hp[1]*v.vector[1]
                res=res0+res1           
            
            #if (res < 1):
            if (res < 1 + self.pertubation1*self.eps*(1/2-random())): # +- k1*eps/2
                v.sign=1 #'-'
            #elif (res > 1+self.eps):
            elif (res > 1+self.eps + self.pertubation1*self.eps*(1/2-random())): # +- k1*eps/2
                v.sign=2 #'+'
            else:
                v.sign=0 #'0'
    
                
    def __set_h_types(self):           
        for h in self.graph.edges:
            h.type=3*h.origin.sign + h.target().sign
            h.twin.type=3*h.target().sign + h.origin.sign
    
    
    def __add_verts(self):
        for h in self.graph.edges[:]:
            if h.type==5 or h.type==7:  #-+ or +-
                if h.type==7:
                    h=h.twin # now h is -+               
                u_vec=h.origin.vector
                w_vec=h.target().vector
                hp=self.matrix[self.iter] 
                hw=hp @ w_vec 
                hd=hp @ (u_vec-w_vec)                      
                #vec = ((1 + self.eps/2 - hw)/ hd) * (u_vec-w_vec) + w_vec
                vec = ((1 + self.eps/2 + self.pertubation2*self.eps*(1/2-random()) - hw)/ hd) * (u_vec-w_vec) + w_vec # +- k2*eps/2
                h1,h2=self.graph.split_edge(h) # h as -+ expected
                h2.origin.vector=vec
                h2.origin.sign=0 # 0
                h1.type=3        # -0
                h1.twin.type=1   # 0-
                h2.type=2        # 0+
                h2.twin.type=6   # +0
                

    def __cut(self):
        """
        cut off "+"-part of graph and draw cycles along the cut
            "+"-part not necessarily connected
        """
        out_cuts=[]
        
        # store all +0 halfedges in queue
        for he in self.graph.edges: 
            if he.type==6: # +0
                if he.face.valid:
                    out_cuts.append(he)
            elif he.type==2: # 0+
                if he.face.valid:
                    out_cuts.append(he.twin)  

                
        # insert edges around V_+ 
        
        # Main idea
        #   stopping points are
        #       case 1: type 2 halfedges
        #       case 2: followers of type 6 halfedges 
        while out_cuts!=[]:
            he=out_cuts[0]
            ho=he.next # go to first stopping point (case 2)
            finished=False
            while True:                
                hi=ho
                while True: 
                    if hi.type==2:
                        self.graph.split_face(hi,ho) # hi and ho are stopping points
                        break
                    hi=hi.next
                # end while: hi is stopping point (case 1)
                ho=hi.next
                while True:
                    if ho.type==6:
                        out_cuts.remove(ho)
                        if ho==he:
                            finished=True
                            break # optional
                        ho=ho.next 
                        self.graph.split_face(ho,hi) # ho and hi are stopping points
                        break
                    else:
                        ho=ho.next
                # end while: ho is stopping point (case 2)
       
                if finished:
                    break
                    
                    
        #a=len(self.graph.get_faces_of_degree(2))          
        
        
        # remove V_+       
        for v in self.graph.vertices[:]:
            if v.sign==2: # +
                while v.halfedge!=None:
                    self.graph.remove_edge(v.halfedge)
                    
        f2=self.graph.get_faces_of_degree(2)
        for f in f2:
            self.graph.remove_edge(f.halfedge)         

                
    def step(self):
        if self.iter < self.matrix.shape[0]:
            t=clock()            
            self.__set_v_signs()            
            self.time1+=clock()-t
            
            t=clock()            
            self.__set_h_types()            
            self.time2+=clock()-t
            
            t=clock()            
            self.__add_verts()            
            self.time3+=clock()-t
            
            t=clock()            
            self.__cut()            
            self.time4+=clock()-t
            
            self.iter+=1
        else:
            print('ApproxVE.step(): no further step to do.')
            
            
    def run(self):
        t=clock()
        while self.iter < self.matrix.shape[0]:
            #print('GraphAlg: Processing inequality {}'.format(self.iter))
            self.step()
        self.time0=clock()-t 
        
        #self.graph.remove_bridges()
        
        self.graph.find_components()
        
    def list_faces(self):
        self.graph.list_faces() 
        
    def list_vertices(self):
        self.graph.list_vertices()
        
    def info(self):
        print('--------------------------------------------------')
        print('Method : GraphAlg')
        print('Problem: {0} with eps= {1}'.format(self.problem.filename,self.eps))
        print('--------------------------------------------------')
        print('CPU-time in seconds .....')
        print('  set_v_signs()  : {}'.format(round(self.time1,3)))
        print('  set_h_types()  : {}'.format(round(self.time2,3)))
        print('  add_verts()    : {}'.format(round(self.time3,3)))
        print('  cut()          : {}'.format(round(self.time4,3)))
        print('  ----------------')
        print('  total          : {}'.format(round(self.time0,3)))
        print('--------------------------------------------------')
        print('Polytope information .....')
        print('  Vertices          : {}'.format(len(self.graph.vertices)))
        print('  Edges             : {}'.format(len(self.graph.edges)))
        print('  Faces             : {}'.format(len(self.graph.faces)))
        print('  Bridges           : {}'.format(len(self.graph.get_bridges())))
        print('  Degree-0 vertices : {}'.format(len(self.graph.get_vertices_of_degree(0))))
        print('  Degree-1 vertices : {}'.format(len(self.graph.get_vertices_of_degree(1))))
        print('  Degree-2 vertices : {}'.format(len(self.graph.get_vertices_of_degree(2))))
        print('  Degree-0 faces    : {}'.format(len(self.graph.get_faces_of_degree(0))))
        print('  Degree-1 faces    : {}'.format(len(self.graph.get_faces_of_degree(1))))
        print('  Degree-2 faces    : {}'.format(len(self.graph.get_faces_of_degree(2))))
        print('--------------------------------------------------')
        nc=len(self.graph.comp_he)
        print('Number of graph components: {}'.format(nc))
        if nc>1:
            for i in range(nc):
                print('  Component {0} has {1} faces'.format(i+1,self.graph.size_of_component(i+1)))
        print('--------------------------------------------------')
        
    def result_to_off(self):
        self.graph.export_to_off()
        
    def kill(self):        
        self.graph.kill()        
        self.graph=None

        

        
