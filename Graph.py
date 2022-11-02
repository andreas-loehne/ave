# -*- coding: utf-8 -*-
"""
Implementation of the half-edge data structure

""" 
from sys import exit
from time import sleep
from random import random,seed
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

H_TYPE={0:'00',1:'0-',2:'0+',3:'-0',4:'--',5:'-+',6:'+0',7:'+-',8:'++'}
V_SIGN={0:'0',1:'-',2:'+'}

class HalfEdge:
    """
    Implements a half-edge
    """
    halfedge_cnt=0
    def __init__(self):
        self.index  = HalfEdge.halfedge_cnt
        HalfEdge.halfedge_cnt+=1
        self.origin = None      #reference to the origin vertex
        self.face   = None      #reference to the incident face
        self.twin   = None      #reference to the twin halfedge
        self.next   = None      #reference to the next halfedge along the boundary of the face
        self.prev   = None      #reference to the preceding halfedge along the boundary of the face
        self.main   = None      #reference to main halfedge in edge (self or twin)
        self.type   = 0         # H-TYPE
        self.comment=''
        
        
    def __str__(self):
        try:
            return 'h({0},{1})'.format(self.origin.index,self.target().index)
        except:
            return 'h(.,.) undefined'
        
        
    def __repr__(self):
        try: 
            return 'h({0},{1})'.format(self.origin.index,self.target().index)
        except:
            return 'h(.,.) undefined.'
        
        
    def follows(self, other):
        self.prev=other
        other.next=self
        
    def followed_by(self, other):
        self.next=other
        other.prev=self
        
    def check_coupling(self):
        p=self.prev
        n=self.next
        t=self.twin
        m=self.main
        if p.next!=self:
            return False
        if n.prev!=self:
            return False
        if t.twin!=self:
            return False
        if m!=self and m!=t:
            return False                
        return True       
        
    def target(self):
        return self.twin.origin        

    def info(self):
        print('HalfEdge')
        print('  Comment:  {}'.format(self.comment))
        print('  Index:    {}'.format(self.index))
        print('  Hash:     {}'.format(self.__hash__()))
        print('  Origin:   {}'.format(self.origin))
        print('  Target:   {}'.format(self.target()))
        print('  Face:     {}'.format(self.face))
        print('  Prev:     {}'.format(self.prev))
        print('  Next:     {}'.format(self.next))
        print('  Main:     {}'.format(self.main))
        print('  Cut:      {}'.format(H_TYPE[self.type]))

        



class Vertex:
    """
    Implements a vertex carrying a vector in R^3 and a sign        
    """
    vertex_cnt=0
    def __init__(self):
        self.index=Vertex.vertex_cnt
        Vertex.vertex_cnt+=1
        self.halfedge=None
        self.vector=None
        self.sign=0 # V_SIGN
        self.face_list=[]
        self.comment=''
        
        
    def __str__(self):
        return 'v({0})'.format(self.index)
        
        
    def __repr__(self):
        return 'v({0})'.format(self.index)
        
        
    def couple(self,he):
        self.halfedge=he
        he.origin=self
        
        
    def __couple_prev(self):
        self.couple(self.halfedge.twin.next)      
        
        
    def __couple_next(self):
        self.couple(self.halfedge.prev.twin)
        
    
    def check_coupling(self):
        if self.halfedge!=None:
            he = self.halfedge
            while True:
                if self.halfedge.origin!= self:
                    return False
                self.__couple_next()
                if self.halfedge==he:
                    break
        return True
    
        
    def update_face_list(self):
        self.face_list=[]
        if self.halfedge!=None:
            he=self.halfedge
            h=he
            while True:                    
                self.face_list.append(h.face)
                h=h.twin.next
                if h==he:
                    break
        
    def get_face_list(self):
        self.update_face_list()
        return self.face_list        
        
    def info(self):
        print('Vertex')
        print('  Comment:  {}'.format(self.comment))
        print('  Index:    {}'.format(self.index))
        print('  Halfedge: {}'.format(self.halfedge))
        print('  Vector:   {}'.format(self.vector))
        print('  Sign:     {}'.format(V_SIGN[self.sign]))     
        
 
 
        
class Face:
    """
    Implements a face         
    """
    face_cnt=0
    def __init__(self):
        self.index=Face.face_cnt
        Face.face_cnt += 1
        self.halfedge=None
        self.component_index=0
        self.vertex_list=[]
        self.comment=''
        self.valid=True
        
        
    def __str__(self):
        return 'f({0})'.format(self.index)
        
        
    def __repr__(self):
        return 'f({0})'.format(self.index)
        
    
    def couple(self,he):
        self.halfedge=he
        he.face=self
        
        
    def __couple_prev(self):
        self.couple(self.halfedge.prev)      
        
        
    def __couple_next(self):
        self.couple(self.halfedge.next)
        
    
    def check_coupling(self):
        he = self.halfedge
        while True:
            if self.halfedge.face!= self:
                return False
            self.__couple_next()
            if self.halfedge==he:
                break
        return True
        

    def update_vertex_list(self):
        self.vertex_list=[]
        if self.halfedge!=None:
            he=self.halfedge
            h=he
            i=0
            while True:  
                i+=1                  
                self.vertex_list.append(h.origin)
                h=h.next
                if h==he:
                    break
            # end while        
        else:
            print('This should not happen') 
        return 0    
            
    def get_vertex_list(self):
        self.update_vertex_list()
        return self.vertex_list
        
        
    def get_neighbors(self,test_index=0,set_index=1):
        nlist=[]
        if self.halfedge==None:
            return nlist
        else:
            he=self.halfedge
            h=he
            while True:
                f=h.twin.face
                if f.component_index==test_index:
                    f.component_index=set_index
                    nlist.append(f)
                h=h.next
                if h==he:
                    break
        return nlist               
        
        
    def info(self):
        print('Face')
        print('  Comment:         {}'.format(self.comment))
        print('  Index:           {}'.format(self.index))
        print('  Halfedge:        {}'.format(self.halfedge))
        print('  component_index: {}:'.format(self.component_index))




class Graph:
    """
    Implements a doubly-connected half-edge list         
    """
    def __init__(self,dim): # initialize as simplex
        self.vertices=[]
        self.faces=[]
        self.edges=[]
        self.comp_he=[] # for storing one halfedge for each component
        if dim==3:
            v0=self.__add_vertex()
            v1=self.__add_vertex()
            v2=self.__add_vertex()
            v3=self.__add_vertex()
            f0=self.__add_face()
            f1=self.__add_face()
            f2=self.__add_face()
            f3=self.__add_face()
            self.__add_edge(v0,v1,f3,f2)
            self.__add_edge(v1,v2,f3,f0)
            self.__add_edge(v2,v0,f3,f1)
            self.__add_edge(v1,v3,f0,f2)
            self.__add_edge(v3,v2,f0,f1)
            self.__add_edge(v0,v3,f2,f1)
            self.__set_prev_next_for_each_halfedge() 
            self.__set_halfedge_for_each_vertex()
            self.__set_halfedge_for_each_face()
        elif dim==2:
            v0=self.__add_vertex()
            v1=self.__add_vertex()
            v2=self.__add_vertex()
            f0=self.__add_face()
            f1=self.__add_face()
            self.__add_edge(v0,v1,f0,f1)
            self.__add_edge(v1,v2,f0,f1)
            self.__add_edge(v2,v0,f0,f1)
            self.__set_prev_next_for_each_halfedge() 
            self.__set_halfedge_for_each_vertex()
            self.__set_halfedge_for_each_face()
            f1.valid=False
        else:
            print('dim must be 2 or 3')
            exit(1)
        self.dim=dim

        
    def __str__(self):
        return 'Graph:\n Halfedges: {0} \n Vertices: {1} \n Faces: {2}'.format(self.edges,self.vertices,self.faces)   
             
        
    def __repr__(self):
        return 'Graph:\n Halfedges: {0} \n Vertices: {1} \n Faces: {2}'.format(self.edges,self.vertices,self.faces) 
               
        
    def __set_prev_next_for_each_halfedge(self):
        """
        slow, for small instances only
        """
        for h0 in self.edges:
            for h in [h0,h0.twin]:
                f1=h.face
                for g0 in self.edges:
                    for g in [g0,g0.twin]:
                        f2=g.face
                        if (f1==f2) and (h.twin.origin == g.origin):
                            g.follows(h)
                            break 
                    
                    
    def __set_halfedge_for_each_vertex(self):
        """
        slow, for small instances only
        """
        for he in self.edges:
            he.origin.halfedge=he
            he.target().halfedge=he.twin
            
            
    def __set_halfedge_for_each_face(self):
        """
        slow, for small instances only
        """
        for he in self.edges:
            he.face.halfedge=he
            he.twin.face.halfedge=he.twin
            
            
    def __add_edge(self,origin,target,face_cw,face_ccw):
        """
        add edge to graph
            .prev and .next to be set manually
        """
        he1=HalfEdge()
        he1.main=he1
        self.edges.append(he1) 
        he2=HalfEdge()
        he2.main=he1
        he1.twin=he2
        he2.twin=he1
        origin.couple(he1)
        target.couple(he2)
        face_cw.couple(he1)
        face_ccw.couple(he2)
        return he1, he2
        
        
    def __merge_faces(self,f1,f2):
        he=f2.halfedge
        f1.halfedge=None # to be set manually
        f2.halfedge=None
        h=he
        while True:
            h.face=f1
            h=h.next
            if h==he:
                break
        f1.valid=f1.valid and f2.valid        
        self.__remove_face(f2)
        
    
        
        
    def __kill_edge(self,he,comment='killed'):
        """
        dirty removal of edge
        """
        self.edges.remove(he.main)
        he.twin.comment=comment
        he.comment=comment
        del he.twin
        del he
        
                
    def __add_vertex(self):
        """
        add vertex to graph
        """
        v=Vertex()
        self.vertices.append(v) 
        return v
        
        
    def __remove_vertex(self,v,comment='removed'):
        """
        remove vertex from graph
        """
        if v.halfedge==None: 
            v.comment=comment     
            self.vertices.remove(v)
            del v   
        else:
            print('Illegal removal of vertex')
            exit(1)    
        
        
    def __add_face(self): 
        """
        add face to graph
        """ 
        f=Face()
        self.faces.append(f) 
        return f
 
        
    def __remove_face(self,f,comment='removed'):
        """
        remove face from graph
        """
        if f.halfedge==None: 
            f.comment=comment     
            self.faces.remove(f)
            del f   
        else:
            print('Illegal removal of face')
            exit(1)
 
        
        
    def split_edge(self,uw):
        """
        clean split edge
        splits uw into edges uv and vw
        """
        wu=uw.twin
        u=uw.origin
        w=wu.origin        
        v=self.__add_vertex() # new vertex
        
        # new edge uv
        uv,vu=self.__add_edge(u,v,uw.face,uw.twin.face) 
        u.couple(uv)

        # edge uw changed to vw 
        vw=uw
        wv=vw.twin
        v.couple(vw)

        if uw.prev==wu:
            uv.follows(vu)
        else:
            tmp1=uw.prev
            tmp2=wu.next
            uv.follows(tmp1)
            vu.followed_by(tmp2)
        uv.followed_by(vw)
        vu.follows(wv)
        return uv,vw


        
    def split_face(self,he1,he2):
        """
        clean splitting of face
        """
        
        if he1.twin!=he2 and he2.next!=he1 and he1.next!=he2: # no parallel edges
        #if True:
            u=he1.origin
            w=he2.origin
            if u!=w:
                f1=he1.face
                f2=self.__add_face()
                uw,wu=self.__add_edge(u,w,f2,f1)
            
                tmp1=he1.prev
                tmp2=he2.prev
                uw.follows(tmp1)
                uw.followed_by(he2)
                wu.follows(tmp2)
                wu.followed_by(he1)
            
                h=he2
                while True:
                    h.face=f2
                    h=h.next
                    if h==uw:
                        break
                        
                        
    def remove_edge(self,he):
        """
        clean removal of edge
        """
        
        f1=he.face
        f2=he.twin.face
        o=he.origin
        t=he.target()
        on=he.twin.next
        op=he.prev
        tn=he.next
        tp=he.twin.prev
        

        if f1!=f2:
            self.__merge_faces(f2,f1)   # f2 survives
            f1=f2                       # just a pointer to f2
        elif tp==he and op==he.twin:
            he.face=None
            he.twin.face=None
            f2.halfedge=None
            self.__remove_face(f2)
        elif op!=he.twin and tp!=he:
            # split graph by inserting a new face
            f1=self.__add_face()
            h=on
            h.face=f1
            while True:
                h=h.next
                h.face=f1
                if h==op:
                    break
        else:
            pass
    

        if tp==he: # t is incident to he.twin only
            t.halfedge=None
            he.twin.origin=None
            self.__remove_vertex(t,'removed in remove_edge 1')
        else:
            tn.follows(tp)
            he.twin.follows(he)
            f2.halfedge=tp
            t.halfedge=tn
        if op==he.twin: # o is incident to he only
            o.halfedge=None
            he.origin=None
            self.__remove_vertex(o,'removed in remove_edge 2')
        else:
            on.follows(op)
            he.follows(he.twin)
            f1.halfedge=op
            o.halfedge=on
        self.__kill_edge(he,'removed')
            
        
    def __update_indices(self):
        """
        updates the vertex indices, so that the maximum index
        matches the number of vertices
        """
        ind = 0
        for vertex in self.vertices:
            vertex.index=ind
            ind+=1  
        
 
    def export_to_off(self,filename='graph.off'):
        """
        write result to OFF-file
        """
        self.__update_indices()
        n_vertices = len(self.vertices)
        #n_faces = len(self.faces)
        n_faces=sum([f.valid for f in self.faces])
        n_edges = len(self.edges)
        file1 = open(filename, 'w')
        file1.write('OFF\n\n')
        file1.write('{0} {1} {2}\n'.format(n_vertices,n_faces,n_edges))
        if self.dim==3:
            for v in self.vertices:
                file1.write("{0} {1} {2}\n".format(v.vector[0],v.vector[1],v.vector[2]))
        else:
            for v in self.vertices:
                file1.write("{0} {1} {2}\n".format(v.vector[0],v.vector[1],0))
        for f in self.faces:
            if f.valid:
                verts=f.get_vertex_list()
                file1.write("{} ".format(len(verts)))
                for v in verts:
                    file1.write("{} ".format(v.index))
                if f.component_index>1:
                    seed(f.component_index)
                    r=random()
                    g=random()
                    b=random()
                    file1.write(" {0} {1} {2} 0.2".format(r,g,b))
                else:
                    file1.write(" 0.4 0.4 0.4 0.2")
                file1.write("\n")
        file1.close()    
        
    def plot(self):        
        fig = plt.gca(projection = '3d')
        fig.axis('off')
        patches = []        
        # Faces
        for f in self.faces:
            if f.valid:
                verts=f.get_vertex_list()
                matrix=[]
                for v in verts:
                    vec=v.vector.tolist()
                    if self.dim==2:
                        vec.append(0)
                    matrix.append(vec)                  
                patches.append(matrix)

        # Edges
        for e in self.edges:
            vec1=e.origin.vector.tolist()
            vec2=e.target().vector.tolist()
            if self.dim==2:
                vec1.append(0)
                vec2.append(0)
            matrix=[vec1,vec2]
            patches.append(matrix)

        if len(patches) >= 1:
            pc = Poly3DCollection(patches,
                                  facecolor = [.9,.9,.9],
                                  edgecolor = [.4,.4,.4],
                                  linestyle = '-',
                                  linewidth=.2,
                                  alpha = 0.9)
            fig.add_collection3d(pc)
            
        # Vertices
        x=[]
        y=[]
        z=[]
        for v in self.vertices:
            x.append(v.vector[0])
            y.append(v.vector[1])
            if self.dim==3:
                z.append(v.vector[2])
            else:
                z.append(0)
            #matrix.append(vec)
            
        fig.scatter(x, y, z,
            s = 0,
            color = [.4,.4,.4],
            marker = 'o',
            linewidth=0)
        
        if self.dim==2:
            fig.view_init(90, -90)  
        plt.show()
        
    def list_faces(self):
        for f in self.faces:
            he=f.halfedge
            h=he
            string='Face {}: '.format(f)
            vl=f.get_vertex_list()
            for v in vl:
                string += str(v)+' '
            print(string)               
            
            
    def list_vertices(self):
        for v in self.vertices:
            he=v.halfedge
            h=he
            string='Vertex {}: '.format(v)
            fl=v.get_face_list()
            for f in fl:
                string += str(f)+' '
            print(string)
            
    
    def find_components(self):
        """
        Computes the components of a graph     
        """
        faces=self.faces.copy()
        i=0
        while faces!=[]:
            i+=1
            f=faces.pop(0)
            self.comp_he.append(f.halfedge) # store one halfedge for each component
            f.component_index=i
            queue=[f]       
            while queue!=[]:
                f=queue.pop(0)
                flist=f.get_neighbors(test_index=0,set_index=i)
                for f in flist:
                    queue.append(f)
                    faces.remove(f)
                    
                    
    def number_of_components(self):
        return(len(self.comp_he))
                           
        
        
    def get_faces_of_component(self,i):
        flist=[]
        for f in self.faces:
            if f.component_index==i:
                flist.append(f)
        return flist
        
        
    def get_edges_of_component(self,i):
        helist=[]
        for he in self.edges:
            if he.face.component_index==i:
                helist.append(he)
        return helist
        
        
    def get_vertices_of_component(self,i):
        vlist=[]
        for v in self.vertices:
            if v.halfedge.face.component_index==i:
                vlist.append(v)
        return vlist
        
        
    def component_to_graph(self,i):
        g=Graph()
        g.vertices=self.get_vertices_of_component(i)
        g.edges=self.get_edges_of_component(i)
        g.faces=self.get_faces_of_component(i)
        return g
        
                
    def size_of_component(self,i):
        j=0
        for f in self.faces:
            if f.component_index==i:
                j+=1
        return j    
        
        
    def get_bridges(self):
        blist=[]
        for he in self.edges:
            if he.face==he.twin.face:
                blist.append(he)
        return blist   
        
        
    def remove_bridges(self):
        blist=self.get_bridges()
        for h in blist[:]:
            self.remove_edge(h)
            
            
    def get_vertices_of_degree(self,deg):
        vlist=[]
        for v in self.vertices:
            if v.halfedge==None:
                if deg==0:
                    vlist.append(v)
            elif deg!=0:                    
                he=v.halfedge
                h=he
                for d in range(deg):
                    h=h.prev.twin
                if h==he:
                    vlist.append(v)
        return vlist 
        
        
    def get_faces_of_degree(self,deg):
        flist=[]
        for f in self.faces:
            if f.halfedge==None:
                if deg==0:
                    flist.append(f)
            elif deg!=0:                    
                he=f.halfedge
                h=he
                for d in range(deg):
                    h=h.next
                if h==he:
                    flist.append(f)
        return flist                   
    
            
    def check_coupling_of_faces(self):
        for f in self.faces:
            if f.check_coupling()==False:
                print('-- corrupt coupling of face {}'.format(f))
                print('\a')
                return False
        print(' * coupling of faces correct')
        return True   
        
    def check_coupling_of_vertices(self):
        for v in self.vertices:
            if v.check_coupling()==False:
                print('-- corrupt coupling of vertex {}'.format(v))
                print('\a')
                return False
        print(' * coupling of vertices correct')
        return True     
        
    def check_coupling_of_halfedges(self):
        for he in self.edges:
            for h in [he,he.twin]:
                if h.check_coupling()==False:
                    print('-- corrupt coupling of halfedge {}'.format(h))
                    print('\a')
                    return False
        print(' * coupling of halfedges correct')
        return True
        
        
    def check_he_links_of_halfedges(self):
        for he in self.edges:
            for h in [he,he.twin]:
                if h.next.main not in self.edges:
                    print('-- invalid link h.next for h= {}'.format(h))
                    print('\a')
                    return False
                if h.prev.main not in self.edges:
                    print('-- invalid link h.prev for h= {}'.format(h))
                    print('\a')
                    return False
                if h.twin.main not in self.edges:
                    print('-- invalid link h.twin for h= {}'.format(h))
                    print('\a')
                    return False                               
        print(' * links of halfedges correct')
        return True
        
        
    def check_face_links_of_halfedges(self):
        face_links=set()
        for he in self.edges:
            for h in [he,he.twin]:
                face_links.add(h.face)
        if set(self.faces)!=face_links:
            print('-- face links of halfedges incorrect')
            print('\a')
            return False
        else:                                    
            print(' * face links of halfedges correct')
            return True
            
            
    def check_vertex_links_of_halfedges(self):
        vertex_links=set()
        for he in self.edges:
            for h in [he,he.twin]:
                vertex_links.add(h.origin)
        if set(self.vertices)!=vertex_links:
            print('-- vertex links of halfedges incorrect')
            print('\a')
            return False
        else:                                    
            print(' * vertex links of halfedges correct')
            return True   
            
            
    def check_euler(self):
        f=len(self.faces)
        e=len(self.edges)
        v=len(self.vertices)
        n=self.number_of_components()
        if v-e+f==2*n:                                  
            print(' * Euler formula (v-e+f=2n) holds')
            return True 
        else:
            print('-- Euler formula (v-e+f=2n) does not hold')
            print('\a')
            return False

        
    def check(self):
        print('\n** checking the graph')
        self.check_coupling_of_faces()
        self.check_coupling_of_vertices()
        self.check_coupling_of_halfedges()
        self.check_he_links_of_halfedges()
        self.check_face_links_of_halfedges()
        self.check_vertex_links_of_halfedges()
        self.check_euler()
        print('**\n')
        
        
    def kill(self):        
        for v in self.vertices[:]:
            v.halfedge=None
            self.__remove_vertex(v)
        Vertex.vertex_cnt=0
        for f in self.faces[:]:
            f.halfedge=None
            f.component_index=0
            self.__remove_face(f)
        Face.face_cnt=0
        for h in self.edges[:]:
            self.__kill_edge(h)