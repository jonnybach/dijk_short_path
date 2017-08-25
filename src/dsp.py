"""
standford online course algorithms - dijkstra's shortest path distances

Notes:
Vertices are labeled as positive integers from 1 to 875714. Every row indicates 
an edge, the vertex label in first column is the tail and the vertex label in 
second column is the head (recall the graph is directed, and the edges are 
directed from the first column vertex to the second column vertex). 
"""

#from merge_sort import sort
from __future__ import print_function
from graphviz import Digraph
from heapq import heappush, heappop, heapify

__author__ = 'jonathan'

k_max_cost = 1000000


class vertex:
    
    def __init__(self, id_num):
        self.edges = dict()
        self.explored = False
        self.id = id_num
        self.to_vtxs = list()
        self.run_cost = k_max_cost
        self.min_cost_edge = 0

    '''
    id property
    '''
    def set_id(self, id_num):
        self.id = id_num

    def get_id(self):
        return self.id

    '''
    leader id property
    '''
    def set_ldr_id(self, ldr_id_num):
        self.leader_id = ldr_id_num
        
    def get_ldr_id(self):
        return self.leader_id   

    '''
    explore property
    '''
    def set_explored(self, has_expld):
        self.explored = has_expld
        
    def get_explored(self):
        return self.explored

    '''
    edge methods
    '''
    def add_edge(self, edge):
        self.edges[edge.get_id()] = edge

    def get_edge(self, edge_num):
        edg = self.edges[edge_num]
        return edg

    def remove_edge(self, edge_num):
        del self.edges[edge_num]

    def get_edges(self):
        return self.edges.values()

    def get_out_edges(self):
        '''
        return dictionary of edges where this vertex is the
        tail
        '''
        out_edgs = dict()
        for edg in self.edges.values():
            if edg.get_tail() is self:
                out_edgs[edg.get_id()] = edg
        return out_edgs

    def get_other_vertex_on_edge(self, edge_num):
        edg = self.get_edge(edge_num)
        vtx = 0
        if (self is edg.v1):
            vtx = edg.v2
        else:
            vtx = edg.v1
        return vtx

    '''
    finish time property
    '''
    def get_finish_time(self):
        return self.finsh_time
    
    def set_finish_time(self, fnsh_time):
        self.finsh_time = fnsh_time

    '''
    edge representation as outgoing arc list
    '''
    def set_arc_heads(self, arc_heads):
        self.to_vtxs.append(arc_heads)
    
    def get_arc_heads(self):    
        return self.to_vtxs
    
    def get_arc_head(self, arc_indx):
        return self.to_vtxs(arc_indx)

    '''
    total running cost property
    ''' 
    def set_cost(self, cost):
        self.run_cost = cost
        
    def get_cost(self):
        return self.run_cost

    '''
    min path components property
    '''
    def set_min_cost_edge(self, edge):
        self.min_cost_edge = edge
        
    def get_min_cost_edge(self):
        return self.min_cost_edge


class edge: 
    
    def __init__(self, vtx1, vtx2, id_num):
        self.v1 = vtx1
        self.v2 = vtx2
        self.id = id_num
        self.cost = 1e12
    
    def set_id(self, id_num):
        self.id = id_num
    
    def get_id(self):
        return self.id
    
    def get_head(self):
        return self.v2

    def get_tail(self):
        return self.v1
    
    def get_v1(self):
        return self.v1
    
    def get_v2(self):
        return self.v2

    def find_vtx(self, vertex_id):
        if (self.v1.get_id() == vertex_id):
            return self.v1
        elif (self.v2.get_id() == vertex_id):
            return self.v2
        else:
            return 0
 
    def swap_vtx(self):
        self.v1, self.v2 = self.v2, self.v1
 
    def set_cost(self, cost):
        self.cost = cost
        
    def get_cost(self):
        return self.cost
 
    def replace_vtx(self, vertex_id, new_vertex):
        if (self.v1.get_id() == vertex_id):
            self.v1 = new_vertex
        elif (self.v2.get_id() == vertex_id):
            self.v2 = new_vertex
        return       

    def print_info(self):
        print("\nedge id: {0:3d}".format(self.id))
        print("  vtx1 id: {0:3d}".format(self.v1.get_id()))
        print("  vtx2 id: {0:3d}".format(self.v2.get_id()))



def print_graph(vertices, edges, graph_name):
    
    grf = Digraph(name=graph_name, comment=graph_name, format="pdf")
    
    for edg in edges.values():
        edg_lbl = "E{0:d}-${1:d}".format(edg.get_id(), edg.get_cost())
        v1_id = "{0:d}".format(edg.get_v1().get_id())
        v2_id = "{0:d}".format(edg.get_v2().get_id())
        grf.edge(v1_id, v2_id, label=edg_lbl)
    
    for v in vertices.values():
        v_id = "{0:d}".format(v.get_id())
        v_lbl = "V{0:d}-${1:d}".format(v.get_id(), v.get_cost())
        grf.node(v_id, label = v_lbl, )

    grf.render(view=True, cleanup=False)


def run_find_short_paths(vertex_in, edges_in, strt_vtx_id, vtxs_to_rprt, vtx_crct_dist):
    '''
    compute shortest paths using dijkstra's greedy criteria on edge length
    '''
    
    vertices, edges = vertex_in, edges_in
    #print_graph(vertices, edges, 'orig_graph')

    #run depth first search but sort vertices by reverse finishing time
    print("--> Computing shortest paths starting from vertex->{0:d} ...".format(strt_vtx_id))

    #initialize source vertex
    source_vtx = vertices[strt_vtx_id]
    source_vtx.set_cost(0)
    source_vtx.set_explored(True)
    
    #create vertex heap based on total costs
    # each item added is a tuple with the (cost, vtx id)
    min_cost_hp = []
    for vtx in vertices.values():
        heappush( min_cost_hp, (vtx.get_cost(), vtx.get_id()) )
    
    #main looping routine, pop out next minimum distance value from the heap
    iter = 0
    while (min_cost_hp):
        min_cost, min_cost_vtx_id = heappop(min_cost_hp)
        min_cost_vtx = vertices[min_cost_vtx_id]
        min_cost_vtx.set_explored(True)
        
        #for each edge coming from this search vertex, compute
        # the true cost to the head vertex from the tail vertex
        for edg in min_cost_vtx.get_out_edges().values():  
            true_cost = min_cost_vtx.get_cost() + edg.get_cost()
            hd_vtx = edg.get_head()
            
            if ( hd_vtx.get_explored() ): continue
            
            #if head vertex's distance cost is larger than true cost (initialized to big number)
            if true_cost <= hd_vtx.get_cost():
                hd_vtx.set_cost(true_cost)
                hd_vtx.set_min_cost_edge(edg)
                
                #update the vertex's cost key in the distance heap
                
                #remove the vertex by id in the heap
                # then swap with last value in list, pop, and heapify
                for indx in range(0, len(min_cost_hp)):
                    vtx_cost, vtx_id = min_cost_hp[indx]
                    if (vtx_id == hd_vtx.get_id()):
                        min_cost_hp[indx] = min_cost_hp[-1]
                        min_cost_hp.pop()
                        heapify(min_cost_hp)
                        break
                
                #add the new vertex based on updated cost
                heappush( min_cost_hp, (hd_vtx.get_cost(), hd_vtx.get_id()) )
            
        #print_graph(vertices, edges, 'w_costs_graph_{0:d}'.format(iter))
        iter += 1

    print("--> Shortest path lengths for each vertex:")

    #print results
    for vtx_id in vtxs_to_rprt:
        vtx = vertices[vtx_id]
        print("\tVtx id {0:d}: {1:d}".format(vtx_id, vtx.get_cost()))

    print("\n--> Answers:")
    for ivtx in range(0,len(vtxs_to_rprt)):
        print("\tVtx id {0:d}: {1:d}".format(vtxs_to_rprt[ivtx], vtx_crct_dist[ivtx]))

def get_graph(filename):
    
    vertices = dict()
    edges = dict()

    print("--> Getting graph from file ...")
    
    #read in numbers from text file
    with open(filename) as f:
        #get all lines of text
        lines = f.readlines()
            
    #instantiate all vertices using only the first row integer in all lines
    for line in lines:
        pieces = line.split()
        this_vtx_id = int(pieces[0])
        this_vtx = vertex(this_vtx_id)
        vertices[this_vtx_id] = this_vtx
            
    #print( "total number of vertices read: {0:d}\n".format(len(vertices)) )
        
    #create all edges
    edge_id = 0 #initialize edge id
    for line in lines:
        pieces = line.split()
            
        #get the main vertex
        this_vtx_id = int(pieces[0])
        this_vtx = vertices[this_vtx_id]
            
        #create edges for each vertex where an edge hasn't been created before
        for pc in pieces[1:]:
            othr_vtx_id, cost = pc.split(',')
            othr_vtx = vertices[int(othr_vtx_id)]
            new_edg = edge(this_vtx, othr_vtx, edge_id)
            new_edg.set_cost(int(cost))
            this_vtx.add_edge(new_edg)
            othr_vtx.add_edge(new_edg)
            edges[edge_id] = new_edg
            edge_id += 1
                
    #print( "\ntotal number of edges created: {0:d}".format(len(edges)) )     
    print("-->     Done reading graph from file ...")
    return (vertices, edges)


if __name__ == "__main__":
    
    vertices = dict()
    edges = dict()
    
    test_case = 0
    
    if test_case == 0:
        vertices, edges = get_graph('dijkstraData.txt')
        strt_vtx_id = 1
        vtx_to_rprt = [7,37,59,82,99,115,133,165,188,197]
        #vtx_to_rprt = [ (vtx) for vtx in vtx_to_rprt]
        vtx_crct_dist = [0] * len(vtx_to_rprt)
        run_find_short_paths(vertices, edges, strt_vtx_id, vtx_to_rprt, vtx_crct_dist)
                
    elif test_case == 1:
        vertices, edges = get_graph('test1.txt')
        strt_vtx_id = 1
        vtx_to_rprt = [1,2,3,4]
        vtx_crct_dist= [0,3,3,5]        
        run_find_short_paths(vertices, edges, strt_vtx_id, vtx_to_rprt, vtx_crct_dist)
        

    elif test_case == 2:
        vertices, edges = get_graph('test2.txt')
        strt_vtx_id = 1
        vtx_to_rprt = [1,2,3,4]        
        vtx_crct_dist= [0,3,4,5]
        run_find_short_paths(vertices, edges, strt_vtx_id, vtx_to_rprt, vtx_crct_dist)
