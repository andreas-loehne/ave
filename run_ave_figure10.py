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
#p.eps=pow(10,0)
#p.eps=pow(10,-1)
#p.eps=pow(10,-2)
p.eps=pow(10,-3)

p.filename='ex/bensolvehedron(3,2).csv'
p.init()
a=GraphAlg(p)
a.run()
a.graph.plot()     
a.kill()












    







            















    








