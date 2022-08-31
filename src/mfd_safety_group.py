#!/usr/bin/env python
# coding: utf-8


import os 
import networkx as nx
import sys
import argparse
import more_itertools


def read_input(graphfile,subpathfile):
    
    graph_data = open(graphfile,'r').read().split('\n')
    subpaths = open(subpathfile,'r').read().split('\n')
    
    i = 0
    listOfGraphs = {}
    listOfSubpaths = {}
    k = 0
    
    while(True):
        if "#" in graph_data[i]:
            i = i+1
            N = int(graph_data[i])
            edges = list()
            while(True):
                i = i+1;
                if "#" in graph_data[i]:
                    break;
                if "" == graph_data[i]:
                    break;
                line = graph_data[i].split(" ")
                edges.append((int(line[0]),int(line[1]),float(line[2])))
                if i >= len(graph_data)-1:
                    break;
            G = {'Nodes':N,'list of edges': edges}
            listOfGraphs[k] = G
            k +=1 
        if i >= len(graph_data)-1:
                    break; 
    i = 0
    k = 0
    
    while(True):
        if i >= len(subpaths):
            break;
        if "#" in subpaths[i]:
            i = i + 1 
            N = int(subpaths[i])
            listOfNodes = [list() for n in range(0,N)]
            for j in range(0,N):
                if(i+j+1 >= len(subpaths)):
                    break
                line = subpaths[i+j+1].split(" ")
                listOfNodes[j] = list(map(int,line))
            i = i + N + 1
        
            subs = {'Amount': N,'list of nodes': listOfNodes}
            listOfSubpaths[k] = subs
            k += 1
    
    for k in range(0,len(listOfGraphs)):
        listOfGraphs[k]['number of testing'] = listOfSubpaths[k]['Amount']
        listOfGraphs[k]['group test'] = listOfSubpaths[k]['list of nodes']
    
    return listOfGraphs;
     

def FD_Algorithm(data):
    
    listOfEdges = data['edges']
    solutionMap = data['graph']
    solutionSet = 0
    Kmin = data['minK']
    solutionWeights = 0

    for i in range(1,len(listOfEdges)+1):
        data = MFD(data,i)
        if data['message'] == "solved":
            solutionSet = data['solution']
            solutionWeights = data['weights']
            break;
    
    return data

def SolveInstances(Graphs):
    
    output = {}
    for g in range(0,len(Graphs)): 
    #for g in range(0,1):    
        f = {}
        Edges = set()
        V = set()
        listOfEdges = Graphs[g]['list of edges']
        for k in range(0,len(listOfEdges)):
            (a,b,c) = (listOfEdges[k])
            Edges.add((a,b))
            V.add(a)
            V.add(b)
            f[a,b] = int(float(c))
            
        
        # creation of graphs
        G = nx.DiGraph()
        G.add_edges_from(Edges,weights = f)
        G.add_nodes_from(V)
        
        # creation of adjacent matrix
        AD_in = {}
        AD_out = {}
        
        for v in V:
            setAdj = set()
            for (i,j) in list(G.out_edges(v)):
                if i != v:
                    setAdj.add(i)
                if j != v:
                    setAdj.add(j)
            
            AD_out[v] = list(setAdj)
            
            setAdj = set()
            for (i,j) in list(G.in_edges(v)):
                if i != v:
                    setAdj.add(i)
                if j != v:
                    setAdj.add(j)
            
            AD_in[v] = list(setAdj)
            
        
        # calculating source, sinks and max flows
        S = [x for x in G.nodes() if G.out_degree(x)>=1 and G.in_degree(x)==0]
        D = [x for x in G.nodes() if G.out_degree(x)==0 and G.in_degree(x)>=1]
        maxW = max(f.values())
        
        
        # definition of data
        
        data = {'edges' : Edges,
                'flows' : f,
                'vertices' : V,
                'graph' : G,
                'Kmax' : len(Edges),
                'weights' : {},
                'sources': S,
                'targets': D,
                'message': {},
                'solution': 0,
                'maxFlow': maxW,
                'adj_in': AD_in,
                'adj_out': AD_out,
                'number of testing': Graphs[g]['number of testing'],
                'group test': Graphs[g]['group test'],
                'minK': 2,
                'runtime': 0,
        }
        
        data = FD_Algorithm(data)        
        
        output[g] = data
        
    return output

def FD_Safety(data,K):
    
    #libraries
    import gurobipy as gp
    from gurobipy import GRB
    
    # calculate the minimal flow decomposition based on such graph
    V = data['vertices']
    E = data['edges']
    f = data['flows']
    W = data['maxFlow']
    S = data['sources']
    D = data['targets']
    N = data['number of testing']
    P = data['group test']
    AD_in = data['adj_in']
    AD_out = data['adj_out']
    
    
    try:
        #create extra sets
        T = [(i,j,k) for (i,j) in E for k in range(0,K)]
        SC = [k for k in range(0,K)]
        R = [r for r in range(0,len(P))]
        
        # Create a new model
        model = gp.Model("MFD")
        model.Params.LogToConsole = 0

        # Create variables
        x = model.addVars(T,vtype=GRB.BINARY, name="x")
        w = model.addVars(SC,vtype=GRB.INTEGER, name="w",lb=0)
        z = model.addVars(T,vtype=GRB.CONTINUOUS, name="z",lb=0)
        r = model.addVars(R,vtype=GRB.BINARY,name="r")
    
        model.setObjective(sum(r[p] for p in range(0,len(P))),GRB.MAXIMIZE,)
       
        # flow conservation
        for k in range(0,K):
            for i in V:
                if i in S:
                    model.addConstr(sum(x[i,j,k] for j in AD_out[i])  == 1) 
                if i in D:
                    model.addConstr(sum(x[j,i,k] for j in AD_in[i]) == 1)
                if i not in S and i not in D:
                    model.addConstr(sum(x[i,j,k] for j in AD_out[i]) - sum(x[j,i,k] for j in AD_in[i]) == 0)

        # flow balance
        for (i,j) in E:
            model.addConstr(f[i,j] == sum(z[i,j,k] for k in range(0,K)))

        # linearization
        for (i,j) in E:
            for k in range(0,K):
                model.addConstr(z[i,j,k] <= W*x[i,j,k])
                model.addConstr(w[k] - (1 - x[i,j,k])*W <= z[i,j,k])
                model.addConstr(z[i,j,k] <= w[k])
                
        
        # group testing constraint
        for p in range(0,len(P)):
            for k in range(0,K):
                edges = list(more_itertools.pairwise(P[p]))
                model.addConstr(sum(x[i,j,k] for (i,j) in edges) <= len(edges) - r[p])

        
        # objective function
        model.optimize()
        
        w_sol = [0]*len(range(0,K))
        x_sol = {}
        paths = [list() for i in range(0,K)]
        gamma = [0]*len(P)
    
        
        if model.status == GRB.OPTIMAL:
            data['message'] = 'solved'
            data['runtime'] = model.Runtime;

            
            for v in model.getVars():
                if 'w' in v.VarName:
                    for k in range(0,K):
                        if str(k) in v.VarName:
                            w_sol[k] = v.x
                
                if 'x' in v.VarName:          
                    for (i,j,k) in T:
                        if str(i)+","+str(j)+","+str(k) in v.VarName:
                            x_sol[i,j,k] = v.x
                            
                if 'r' in v.VarName:
                    for r in P:
                        if str(r) in v.VarName:
                            gamma[r] = v.x
                
            for(i,j,k) in T:
                if x_sol[i,j,k] == 1:
                    paths[k].append((i,j))
                    
            
        if model.status == GRB.INFEASIBLE:
            data['message'] = 'unsolved'
        
        data['gamma'] = gamma
        
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

    except AttributeError:
        print('Encountered an attribute error')
    
    return data;

def MFD(data,K):
    
    #libraries
    import gurobipy as gp
    from gurobipy import GRB
    
    # calculate the minimal flow decomposition based on such graph
    V = data['vertices']
    E = data['edges']
    f = data['flows']
    W = data['maxFlow']
    S = data['sources']
    D = data['targets']
    AD_in = data['adj_in']
    AD_out = data['adj_out']
    
    
    try:
        #create extra sets
        T = [(i,j,k) for (i,j) in E for k in range(0,K)]
        SC = [k for k in range(0,K)]
       
        
        # Create a new model
        model = gp.Model("MFD")
        model.Params.LogToConsole = 0

        # Create variables
        x = model.addVars(T,vtype=GRB.BINARY, name="x")
        w = model.addVars(SC,vtype=GRB.INTEGER, name="w",lb=0)
        z = model.addVars(T,vtype=GRB.CONTINUOUS, name="z",lb=0)
    
        model.setObjective(GRB.MINIMIZE)
       
        # flow conservation
        for k in range(0,K):
            for i in V:
                if i in S:
                    model.addConstr(sum(x[i,j,k] for j in AD_out[i])  == 1) 
                if i in D:
                    model.addConstr(sum(x[j,i,k] for j in AD_in[i]) == 1)
                if i not in S and i not in D:
                    model.addConstr(sum(x[i,j,k] for j in AD_out[i]) - sum(x[j,i,k] for j in AD_in[i]) == 0)

        # flow balance
        for (i,j) in E:
            model.addConstr(f[i,j] == sum(z[i,j,k] for k in range(0,K)))

        # linearization
        for (i,j) in E:
            for k in range(0,K):
                model.addConstr(z[i,j,k] <= W*x[i,j,k])
                model.addConstr(w[k] - (1 - x[i,j,k])*W <= z[i,j,k])
                model.addConstr(z[i,j,k] <= w[k])
                    
        
        # objective function
        model.optimize()
        w_sol = [0]*len(range(0,K))
        x_sol = {}
        paths = [list() for i in range(0,K)]
    
        
        if model.status == GRB.OPTIMAL:
            data['message'] = 'solved'
            data['runtime'] = model.Runtime;

            for v in model.getVars():
                if 'w' in v.VarName:
                    for k in range(0,K):
                        if str(k) in v.VarName:
                            w_sol[k] = v.x
                
                if 'x' in v.VarName:          
                    for (i,j,k) in T:
                        if str(i)+","+str(j)+","+str(k) in v.VarName:
                            x_sol[i,j,k] = v.x
                
            for(i,j,k) in T:
                if x_sol[i,j,k] == 1:
                    paths[k].append((i,j))
            
            for k in range(0,len(paths)):
                paths[k] = sorted(paths[k])
                
            data['weights'] = w_sol
            data['solution'] = paths
        
        if model.status == GRB.INFEASIBLE:
            data['message'] = 'unsolved'
            
        
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

    except AttributeError:
        print('Encountered an attribute error')
    
    return data;


def SolveInstanceSafety(data,fp):

    mFDs = SolveInstances(data)   
    
    output = open(fp+".res",'w+')
    
    for g in range(0,len(mFDs)):
            
        print("# graph " + str(g),file=output)
        data = mFDs[g]
        paths = data['solution']
        K = len(data['weights'])
        
        data = FD_Safety(data,K)
        
        # print path
        notSafePath = data['gamma']
        Paths = data['group test']
        for k in range(len(notSafePath)):
            if notSafePath[k] == 1:
                print(k,*Paths[k],file=output)

    return 0;


print("Start")

parser = argparse.ArgumentParser(
    description="""
    Decompose a network flow into a minimum number of weighted paths. 
    This script uses the Gurobi ILP solver.
    """,
    formatter_class=argparse.RawTextHelpFormatter
    )
parser.add_argument('-wt', '--weighttype', type=str, default='int+', help='Type of path weights (default int+):\n   int+ (positive non-zero ints), \n   float+ (positive non-zero floats).')
parser.add_argument('-t', '--threads', type=int, default=0, help='Number of threads to use for the Gurobi solver; use 0 for all threads (default 0).')

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-i', '--input', type=str, help='Input filename', required=True)
requiredNamed.add_argument('-s', '--subpath', type=str, help='Output filename', required=True)
requiredNamed.add_argument('-o', '--output', type=str, help='Output filename', required=True)


args = parser.parse_args()

threads = args.threads
if threads == 0:
    threads = os.cpu_count()
print(f"INFO: Using {threads} threads for the Gurobi solver")

data = read_input(args.input,args.subpath)
data = SolveInstanceSafety(data,args.output)


print("Done")

