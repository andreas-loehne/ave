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

filename='./fig_ave/ave_fig_6b.csv'

file1 = open(filename, 'w')
file1.write("eps,time_addm,verts_addm\n")


for j in range(15):
    print(j)
    p=Problem()
    p.eps=pow(10,-j-1/2)
    p.filename='ex/simplex(3)_recursion_P=P+polar(P)_4.csv'    
    p.init()    
    b=ApproxDDM(p)
    b.run()
    file1.write("{0},{1},{2}\n".format(p.eps,b.time0,len(b.nodes)))  
    b.kill()

    
file1.close()












    







            















    








