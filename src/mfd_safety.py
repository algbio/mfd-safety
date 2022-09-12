#!/usr/bin/env python
# coding: utf-8

import os
import sys
import argparse
import networkx as nx
import gurobipy as gp
from gurobipy import GRB


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

    data['message'] == "unsolved"
    for i in range(data['min_k'], data['max_k'] + 1):
        data = fd_fixed_size(data, i)
        if data['message'] == "solved":
            return data

    return data


def solve_instances(graphs):

    output = list()
    for graph in graphs:

        # creation of NetworkX Graph
        ngraph = nx.DiGraph()
        ngraph.add_weighted_edges_from(graph['edges'])
        vertices = ngraph.nodes
        edges = ngraph.edges
        weights = {
            (u, v): ngraph.edges[u, v]['weight'] for u, v in ngraph.edges
        }

        # creation of adjacency lists
        adj_in = {
            v: set(ngraph.pred[v].keys()) for v in vertices
        }
        adj_out = {
            v: set(ngraph.succ[v].keys()) for v in vertices
        }

        # calculating source, sinks
        sources = [x for x in vertices if ngraph.in_degree(x) == 0]
        sinks = [x for x in vertices if ngraph.out_degree(x) == 0]

        # definition of data
        data = {
            'edges': edges,
            'flows': weights,
            'vertices': vertices,
            'graph': ngraph,
            'max_k': len(edges),
            'weights': list(),
            'sources': sources,
            'sinks': sinks,
            'message': 'unsolved',
            'solution': list(),
            'max_flow_value': max(weights.values()),
            'adj_in': adj_in,
            'adj_out': adj_out,
            'min_k': 1,
            'runtime': 0,
        }

        output.append(mfd_algorithm(data))

    return output


def build_base_ilp_model(data, size):

    vertices = data['vertices']
    edges = data['edges']
    flow = data['flows']
    max_flow_value = data['max_flow_value']
    sources = data['sources']
    sinks = data['sinks']
    adj_in = data['adj_in']
    adj_out = data['adj_out']

    # create extra sets
    T = [(u, v, k) for (u, v) in edges for k in range(size)]
    SC = list(range(size))

    # Create a new model
    model = gp.Model("MFD")
    model.Params.LogToConsole = 0

    # Create variables
    x = model.addVars(T, vtype=GRB.BINARY, name="x")
    w = model.addVars(SC, vtype=GRB.INTEGER, name="w", lb=0)
    z = model.addVars(T, vtype=GRB.CONTINUOUS, name="z", lb=0)

    model.setObjective(GRB.MINIMIZE)

    # flow conservation
    for k in range(size):
        for v in vertices:
            if v in sources:
                model.addConstr(sum(x[v, w, k] for w in adj_out[v]) == 1)
            if v in sinks:
                model.addConstr(sum(x[u, v, k] for u in adj_in[v]) == 1)
            if v not in sources and v not in sinks:
                model.addConstr(sum(x[v, w, k] for w in adj_out[v]) - sum(x[u, v, k] for u in adj_in[v]) == 0)

    # flow balance
    for (u, v) in edges:
        model.addConstr(flow[u, v] == sum(z[u, v, k] for k in range(size)))

    # linearization
    for (u, v) in edges:
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

    if model.status == GRB.OPTIMAL:
        edges = data['edges']
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


def compute_in_out_flows(mfds):
    for mfd in mfds:
        vertices = list(mfd['graph'].nodes)
        mfd['in_flow'] = {v: sum([mfd['flows'][u, v] for u in mfd['adj_in'][v]]) for v in vertices}
        mfd['out_flow'] = {u: sum([mfd['flows'][u, v] for v in mfd['adj_out'][u]]) for u in vertices}

    return mfds


def solve_instances_safety(graphs, output_file):

    mfds = solve_instances(graphs)
    output = open(output_file, 'w+')
    output_counters = open(f"{output_file}.count", 'w+')

    for g, mfd in enumerate(mfds):

        output.write(f"# graph {g}\n")
        output_counters.write(f"# graph {g}\n")
        paths = mfd['solution']

        ilp_count = 0  # Number of ILP calls per graph

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
                        maximal_safe_paths.append(path[first:])
                    break

                ilp_count += 1
                if is_safe(mfd, len(paths), path[first:last + 2]):
                    last = last + 1
                    extending = True
                else:
                    if extending:
                        maximal_safe_paths.append(path[first:last + 1])
                        extending = False
                    first = first + 1

            # print path
            for max_safe in maximal_safe_paths:
                output.write('-1 ')
                output.write(' '.join(map(str, sorted(set([item for t in max_safe for item in t])))))
                output.write('\n')

            output_counters.write(f'{ilp_count}\n')

    output.close()
    output_counters.close()


def solve_instances_safety_using_excess_flow(graphs, output_file):

    mfds = solve_instances(graphs)
    mfds = compute_in_out_flows(mfds)
    output = open(output_file, 'w+')
    output_counters = open(f"{output_file}.count", 'w+')

    for g, mfd in enumerate(mfds):

        output.write(f"# graph {g}\n")
        output_counters.write(f"# graph {g}\n")
        paths = mfd['solution']

        ilp_count = 0  # Number of ILP calls per graph

        for path in paths:

            maximal_safe_paths = list()
            excess_flow = mfd['flows'][path[0]]

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
                        maximal_safe_paths.append(path[first:])
                    break

                extended_excess_flow = excess_flow - mfd['out_flow'][path[last][1]] + mfd['flows'][path[last+1]]
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
                        maximal_safe_paths.append(path[first:last + 1])
                        extending = False

                    excess_flow += mfd['in_flow'][path[first][1]] - mfd['flows'][path[first]]
                    first = first + 1

            # print path
            for max_safe in maximal_safe_paths:
                output.write('-1 ')
                output.write(' '.join(map(str, sorted(set([item for t in max_safe for item in t])))))
                output.write('\n')

            output_counters.write(f'{ilp_count}\n')

    output.close()
    output_counters.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="""
        Computes maximal safe paths for Minimum Flow Decomposition. 
        This script uses the Gurobi ILP solver.
        """,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-wt', '--weighttype', type=str, default='int+',
                        help='Type of path weights (default int+):\n   int+ (positive non-zero ints), \n   float+ (positive non-zero floats).')
    parser.add_argument('-t', '--threads', type=int, default=0,
                        help='Number of threads to use for the Gurobi solver; use 0 for all threads (default 0).')
    parser.add_argument('-uef', '--use-excess-flow', action='store_true', help='Use excess flow of a path to save ILP calls')

    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-i', '--input', type=str, help='Input filename', required=True)
    requiredNamed.add_argument('-o', '--output', type=str, help='Output filename', required=True)

    args = parser.parse_args()

    threads = args.threads
    if threads == 0:
        threads = os.cpu_count()
    print(f"INFO: Using {threads} threads for the Gurobi solver")

    if args.use_excess_flow:
        solve_instances_safety_using_excess_flow(read_input(args.input), args.output)
    else:
        solve_instances_safety(read_input(args.input), args.output)
