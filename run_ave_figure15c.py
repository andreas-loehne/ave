#import ApproxVE
import numpy as np
import Graph
import Div
from random import seed
from random import random
from time import sleep
from Problem import Problem
from GraphAlg import GraphAlg
from ApproxDDM import ApproxDDM

filename='./fig_ave/ave_fig_5.csv'

file1 = open(filename, 'w')
file1.write("time_ga,time_addm,verts_ga,verts_addm\n")


for j in range(6):
    p=Problem()
    p.eps=1e-8
    p.filename='ex/simplex(3)_recursion_P=P+polar(P)_{0}.csv'.format(j+1)    
    p.init()
    
    a=GraphAlg(p)
    a.run()
    
    a.info()
    
    b=ApproxDDM(p)
    b.run()
     
    file1.write("{0},{1},{2},{3}\n".format(a.time0,b.time0,len(a.graph.faces)+len(a.graph.vertices),len(b.facets)+len(b.nodes)))  

    a.kill()
    b.kill()
    
file1.close()












    







            















    








