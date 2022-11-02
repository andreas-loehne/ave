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



p=Problem()
p.eps=pow(10,-9)

#p.filename='ex/simplex(3)_recursion_P=P+polar(P)_1.csv'
#p.filename='ex/simplex(3)_recursion_P=P+polar(P)_2.csv'
#p.filename='ex/simplex(3)_recursion_P=P+polar(P)_3.csv'
#p.filename='ex/simplex(3)_recursion_P=P+polar(P)_4.csv'
#p.filename='ex/simplex(3)_recursion_P=P+polar(P)_5.csv'
p.filename='ex/simplex(3)_recursion_P=P+polar(P)_6.csv'
p.init()
a=GraphAlg(p)
a.run()
a.graph.plot()     
a.kill()












    







            















    








