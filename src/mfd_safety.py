#!/usr/bin/env python
# coding: utf-8

import os
import sys
import argparse
import networkx as nx
import gurobipy as gp
from gurobipy import GRB
from collections import deque


def get_edge(raw_edge):

    parts = raw_edge.split()
    return int(parts[0]), int(parts[1]), float(parts[2])


def get_graph(raw_graph):

    graph = {
        'n': 0,
        'edges': list()
    }

    try:
        lines = raw_graph.split('\n')[1:]
        if not lines[-1]:
            lines = lines[:-1]
        graph['n'], graph['edges'] = int(lines[0]), [get_edge(raw_e) for raw_e in lines[1:]]

    finally:
        return graph


def read_input_graphs(graph_file):

    graphs_raw = open(graph_file, 'r').read().split('#')[1:]
    return [get_graph(raw_g) for raw_g in graphs_raw]


def read_input(graph_file):

    return read_input_graphs(graph_file)


def mfd_algorithm(data):

    data['message'] = 'unsolved'
    for i in range(1, len(data['graph'].edges) + 1):
        data = fd_fixed_size(data, i)
        if data['message'] == 'solved':
            return data

    return data


def solve_instances(graphs_metadata):

    return [mfd_algorithm(graph_metadata) for graph_metadata in graphs_metadata]


def build_base_ilp_model(data, size):

    graph = data['graph']
    max_flow_value = data['max_flow_value']
    sources = data['sources']
    sinks = data['sinks']

    # create extra sets
    T = [(u, v, k) for (u, v) in graph.edges for k in range(size)]
    SC = list(range(size))

    # Create a new model
    model = gp.Model('MFD')
    model.Params.LogToConsole = 0

    # Create variables
    x = model.addVars(T, vtype=GRB.BINARY, name='x')
    w = model.addVars(SC, vtype=GRB.INTEGER, name='w', lb=0)
    z = model.addVars(T, vtype=GRB.CONTINUOUS, name='z', lb=0)

    model.setObjective(GRB.MINIMIZE)

    # flow conservation
    for k in range(size):
        for v in graph.nodes:
            if v in sources:
                model.addConstr(sum(x[v, w, k] for w in graph.successors(v)) == 1)
            if v in sinks:
                model.addConstr(sum(x[u, v, k] for u in graph.predecessors(v)) == 1)
            if v not in sources and v not in sinks:
                model.addConstr(sum(x[v, w, k] for w in graph.successors(v)) - sum(x[u, v, k] for u in graph.predecessors(v)) == 0)

    # flow balance
    for (u, v, f) in graph.edges(data='flow'):
        model.addConstr(f == sum(z[u, v, k] for k in range(size)))

    # linearization
    for (u, v) in graph.edges:
        for k in range(size):
            model.addConstr(z[u, v, k] <= max_flow_value * x[u, v, k])
            model.addConstr(w[k] - (1 - x[u, v, k]) * max_flow_value <= z[u, v, k])
            model.addConstr(z[u, v, k] <= w[k])

    return model, x, w, z


def build_ilp_model_avoiding_path(data, size, path):

    model, x, _, _ = build_base_ilp_model(data, size)

    # safety test constraint
    for k in range(size):
        model.addConstr(sum(x[u, v, k] for (u, v) in path) <= len(path) - 1)

    return model


def get_solution(model, data, size):

    data['weights'], data['solution'] = list(), list()

    if model.status == GRB.OPTIMAL:
        edges = data['graph'].edges
        T = [(u, v, k) for (u, v) in edges for k in range(size)]

        w_sol = [0] * len(range(size))
        paths = [list() for _ in range(size)]
        for k in range(size):
            w_sol[k] = model.getVarByName(f'w[{k}]').x
        for (u, v, k) in T:
            if model.getVarByName(f'x[{u},{v},{k}]').x == 1:
                paths[k].append((u, v))
        for k in range(len(paths)):
            paths[k] = sorted(paths[k])

        data['weights'], data['solution'] = w_sol, paths

    return data


def update_status(data, model):

    if model.status == GRB.OPTIMAL:
        data['message'] = 'solved'
        data['runtime'] = model.Runtime

    if model.status == GRB.INFEASIBLE:
        data['message'] = 'unsolved'
        data['runtime'] = 0

    return data


def fd_fixed_size_forbidding_path(data, size, path):

    # calculate a flow decomposition into size paths avoiding path
    try:
        # Create a new model
        model = build_ilp_model_avoiding_path(data, size, path)
        model.optimize()
        data = update_status(data, model)

    except gp.GurobiError as e:
        print(f'Error code {e.errno}: {e}', file=sys.stderr)

    except AttributeError:
        print('Encountered an attribute error', file=sys.stderr)

    return data


def is_safe(mfd, k, path):
    return fd_fixed_size_forbidding_path(mfd, k, path)['message'] != 'solved'


def fd_fixed_size(data, size):

    # calculate a flow decomposition into size paths
    try:
        # Create a new model
        model, _, _, _ = build_base_ilp_model(data, size)

        # objective function
        model.optimize()
        data = update_status(data, model)
        data = get_solution(model, data, size)

    except gp.GurobiError as e:
        print(f'Error code {e.errno}: {e}', file=sys.stderr)

    except AttributeError:
        print('Encountered an attribute error', file=sys.stderr)

    return data


def output_maximal_safe_path(output, max_safe):
    output.write('-1 ')
    output.write(' '.join(map(str, sorted(set([item for t in max_safe for item in t])))))
    output.write('\n')


def compute_maximal_safe_paths(mfd):

    ilp_count = 0  # Number of ILP calls per graph
    paths = mfd['solution']
    max_safe_paths = list()

    for path in paths:

        maximal_safe_paths = list()

        # two-finger algorithm
        # Invariant: path[first:last+1] (python notation) is safe
        first = 0
        last = 0
        extending = True

        while True:

            if first > last:
                last = first
                extending = True

            if last >= len(path) - 1:
                if extending:
                    maximal_safe_paths.append((first, len(path)))
                break

            ilp_count += 1
            if is_safe(mfd, len(paths), path[first:last + 2]):
                last = last + 1
                extending = True
            else:
                if extending:
                    maximal_safe_paths.append((first, last + 1))
                    extending = False
                first = first + 1

        max_safe_paths.append(maximal_safe_paths)

    return paths, max_safe_paths, ilp_count


def get_flow(e, mfd):
    return mfd['graph'].edges[e]['flow']


def compute_maximal_safe_paths_using_excess_flow(mfd):

    ilp_count = 0  # Number of ILP calls per graph
    paths = mfd['solution']
    max_safe_paths = list()

    for path in paths:

        maximal_safe_paths = list()
        excess_flow = get_flow(path[0], mfd)

        # two-finger algorithm
        # Invariant: path[first:last+1] (python notation) is safe
        first = 0
        last = 0
        extending = True

        while True:

            if first > last:
                last = first
                extending = True

            if last >= len(path) - 1:
                if extending:
                    maximal_safe_paths.append((first, len(path)))
                break

            extended_excess_flow = excess_flow - mfd['out_flow'][path[last][1]] + get_flow(path[last + 1], mfd)
            if extended_excess_flow > 0:
                last = last + 1
                extending = True
                excess_flow = extended_excess_flow
            elif is_safe(mfd, len(paths), path[first:last + 2]):
                ilp_count += 1
                last = last + 1
                extending = True
                excess_flow = extended_excess_flow
            else:
                ilp_count += 1
                if extending:
                    maximal_safe_paths.append((first, last + 1))
                    extending = False

                excess_flow += mfd['in_flow'][path[first][1]] - get_flow(path[first], mfd)
                first = first + 1

        max_safe_paths.append(maximal_safe_paths)

    return paths, max_safe_paths, ilp_count


def get_out_tree(ngraph, v, processed):

    # Get root of out_tree
    root = v
    while ngraph.in_degree(root) == 1:
        root = next(ngraph.predecessors(root))

    leaf_edges = list()
    edges = deque(ngraph.out_edges(root))


    while edges:
        u, v = edges.popleft()
        processed[u] = True
        if ngraph.in_degree(v) != 1:
            leaf_edges.append((u, v))
        else:
            for e in ngraph.out_edges(v):
                edges.append(e)

    # Filter out uncompressible edges
    leaf_edges = [(u, v) for (u, v) in leaf_edges if u != root]

    return root, leaf_edges


def get_out_trees(ngraph):

    out_trees = dict()
    processed = {v: ngraph.in_degree(v) != 1 for v in ngraph.nodes}
    for v in ngraph.nodes:
        if not processed[v]:
            root, leaf_edges = get_out_tree(ngraph, v, processed)
            out_trees[root] = leaf_edges
    return out_trees


def get_y_to_v_contraction(ngraph):

    out_trees = get_out_trees(ngraph)

    mngraph = nx.MultiDiGraph()
    for u, v in ngraph.edges():
        if ngraph.in_degree(v) != 1 and ngraph.in_degree(u) != 1:
            mngraph.add_edge(u, v, flow=ngraph.edges[u, v]['flow'], pred=u)

    for root, leaf_edges in out_trees.items():
        for u, v in leaf_edges:
            mngraph.add_edge(root, v, flow=ngraph.edges[u, v]['flow'], pred=u)

    return mngraph


def get_expanded_path(path, mfd):
    return path


def compute_graphs_metadata(graphs, use_excess_flow, use_y_to_v):

    output = list()
    for graph in graphs:

        # creation of NetworkX Graph
        ngraph = nx.DiGraph()
        ngraph.add_weighted_edges_from(graph['edges'], weight='flow')
        if use_y_to_v:
            original = ngraph
            ngraph = get_y_to_v_contraction(ngraph)
            ngraph = original  # Still not working with the multigraph

        # calculating source, sinks
        sources = [x for x in ngraph.nodes if ngraph.in_degree(x) == 0]
        sinks = [x for x in ngraph.nodes if ngraph.out_degree(x) == 0]

        if use_excess_flow:
            in_flow = {v: sum([f for u, v, f in ngraph.in_edges(v, data='flow')]) for v in ngraph.nodes}
            out_flow = {u: sum([f for u, v, f in ngraph.out_edges(u, data='flow')]) for u in ngraph.nodes}

        # definition of data
        data = {
            'graph': ngraph,
            'sources': sources,
            'sinks': sinks,
            'max_flow_value': max(ngraph.edges(data='flow'), key=lambda e: e[-1])[-1],
            'compute': compute_maximal_safe_paths_using_excess_flow if use_excess_flow else compute_maximal_safe_paths,
            'in_flow': in_flow if use_excess_flow else None,
            'out_flow': out_flow if use_excess_flow else None,
            'original': original if use_y_to_v else None
        }

        output.append(data)

    return output


def solve_instances_safety(graphs, output_file, use_excess_flow, use_y_to_v):

    graphs = compute_graphs_metadata(graphs, use_excess_flow, use_y_to_v)
    mfds = solve_instances(graphs)

    output = open(output_file, 'w+')
    output_counters = open(f'{output_file}.count', 'w+')

    for g, mfd in enumerate(mfds):

        output.write(f'# graph {g}\n')
        output_counters.write(f'# graph {g}\n')

        paths, max_safe_paths, ilp_count = mfd['compute'](mfd)

        for path, maximal_safe_paths in zip(paths, max_safe_paths):
            # print path
            for i, j in maximal_safe_paths:
                max_path = path[i:j]
                if use_y_to_v:
                    max_path = get_expanded_path(max_path, mfd)
                output_maximal_safe_path(output, max_path)

        output_counters.write(f'{ilp_count}\n')

    output.close()
    output_counters.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='''
        Computes maximal safe paths for Minimum Flow Decomposition. 
        This script uses the Gurobi ILP solver.
        ''',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-wt', '--weighttype', type=str, default='int+',
                        help='Type of path weights (default int+):\n   int+ (positive non-zero ints), \n   float+ (positive non-zero floats).')
    parser.add_argument('-t', '--threads', type=int, default=0,
                        help='Number of threads to use for the Gurobi solver; use 0 for all threads (default 0).')
    parser.add_argument('-uef', '--use-excess-flow', action='store_true', help='Use excess flow of a path to save ILP calls')
    parser.add_argument('-uy2v', '--use-y-to-v', action='store_true', help='Use Y to V contraction of the input graphs')

    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-i', '--input', type=str, help='Input filename', required=True)
    requiredNamed.add_argument('-o', '--output', type=str, help='Output filename', required=True)

    args = parser.parse_args()

    threads = args.threads
    if threads == 0:
        threads = os.cpu_count()
    print(f'INFO: Using {threads} threads for the Gurobi solver')

    solve_instances_safety(read_input(args.input), args.output, args.use_excess_flow, args.use_y_to_v)
