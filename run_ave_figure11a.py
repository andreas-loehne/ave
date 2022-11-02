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

filename='./fig_ave/ave_fig_3.csv'

file1 = open(filename, 'w')
file1.write("eps,time_ga,time_addm\n")


for j in range(14):
    p=Problem()
    p.eps=pow(10,-j/4)
    p.filename='ex/bensolvehedron(3,2).csv'    
    p.init()
    
    a=GraphAlg(p)
    a.run()
    
    b=ApproxDDM(p)
    b.run()
     
    file1.write("{0},{1},{2}\n".format(p.eps,a.time0,b.time0))  
    
    Div.SubsetTest(a.graph.vertices,b.nodes)

    a.kill()
    b.kill()
    
file1.close()












    







            















    








