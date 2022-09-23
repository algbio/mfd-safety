#!/usr/bin/env python
# coding: utf-8

import os
import sys
import argparse
import networkx as nx
import gurobipy as gp
from gurobipy import GRB
from collections import deque
from bisect import bisect

ilp_counter = 0


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


def build_base_ilp_model(data, size):

    graph = data['graph']
    max_flow_value = data['max_flow_value']
    sources = data['sources']
    sinks = data['sinks']

    # create extra sets
    T = [(u, v, i, k) for (u, v, i) in graph.edges(keys=True) for k in range(size)]
    SC = list(range(size))

    # Create a new model
    model = gp.Model('MFD')
    model.Params.LogToConsole = 0

    # Create variables
    x = model.addVars(T, vtype=GRB.BINARY, name='x')
    w = model.addVars(SC, vtype=GRB.INTEGER, name='w', lb=0)
    z = model.addVars(T, vtype=GRB.CONTINUOUS, name='z', lb=0)

    # flow conservation
    for k in range(size):
        for v in graph.nodes:
            if v in sources:
                model.addConstr(sum(x[v, w, i, k] for _, w, i in graph.out_edges(v, keys=True)) == 1)
            if v in sinks:
                model.addConstr(sum(x[u, v, i, k] for u, _, i in graph.in_edges(v, keys=True)) == 1)
            if v not in sources and v not in sinks:
                model.addConstr(sum(x[v, w, i, k] for _, w, i in graph.out_edges(v, keys=True)) - sum(x[u, v, i, k] for u, _, i in graph.in_edges(v, keys=True)) == 0)

    # flow balance
    for (u, v, i, f) in graph.edges(keys=True, data='flow'):
        model.addConstr(f == sum(z[u, v, i, k] for k in range(size)))

    # linearization
    for (u, v, i) in graph.edges(keys=True):
        for k in range(size):
            model.addConstr(z[u, v, i, k] <= max_flow_value * x[u, v, i, k])
            model.addConstr(w[k] - (1 - x[u, v, i, k]) * max_flow_value <= z[u, v, i, k])
            model.addConstr(z[u, v, i, k] <= w[k])

    return model, x, w, z


def build_ilp_model_avoiding_path(data, size, path):

    model, x, _, _ = build_base_ilp_model(data, size)

    # safety test constraint
    for k in range(size):
        model.addConstr(sum(x[u, v, i, k] for (u, v, i) in path) <= len(path) - 1)

    return model


def get_solution(model, data, size):

    data['weights'], data['solution'] = list(), list()

    if model.status == GRB.OPTIMAL:
        graph = data['graph']
        T = [(u, v, i, k) for (u, v, i) in graph.edges(keys=True) for k in range(size)]

        w_sol = [0] * len(range(size))
        paths = [list() for _ in range(size)]
        for k in range(size):
            w_sol[k] = model.getVarByName(f'w[{k}]').x
        for (u, v, i, k) in T:
            if model.getVarByName(f'x[{u},{v},{i},{k}]').x == 1:
                paths[k].append((u, v, i))
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

        global ilp_counter
        ilp_counter += 1

        data = update_status(data, model)

    except gp.GurobiError as e:
        print(f'Error code {e.errno}: {e}', file=sys.stderr)

    except AttributeError:
        print('Encountered an attribute error', file=sys.stderr)

    return data


def is_safe(mfd, k, path, i, j, safe_database=None, unsafe_database=None):
    '''
    --> safe_database[i] stores the safe_interval in path such that its left endpoint is the predecessor of i in the database.
    We assume that safe_database contains a maximal set of intervals.
    --> unsafe_database[j] stores the unsafe_interval in path such that its right endpoint is the predecesor of j in the database.
    We assume that unsafe_database contains a minimal set of intervals.
    '''

    if safe_database:
        si, sj = safe_database[i]
        if si <= i and j <= sj:
            return True
    if unsafe_database:
        ui, uj = unsafe_database[j]
        if i <= ui and uj <= j:
            return False

    return fd_fixed_size_forbidding_path(mfd, k, path[i:j])['message'] != 'solved'


def fd_fixed_size(data, size):

    # calculate a flow decomposition into size paths
    try:
        # Create a new model
        model, _, _, _ = build_base_ilp_model(data, size)

        # objective function
        model.optimize()

        global ilp_counter
        ilp_counter += 1

        data = update_status(data, model)
        data = get_solution(model, data, size)

    except gp.GurobiError as e:
        print(f'Error code {e.errno}: {e}', file=sys.stderr)

    except AttributeError:
        print('Encountered an attribute error', file=sys.stderr)

    return data


def output_maximal_safe_path(output, max_safe):
    output.write('-1 ')
    output.write(' '.join(map(str, sorted(set([item for t in max_safe for item in t[:-1]])))))
    output.write('\n')


def get_flow(e, mfd):

    return mfd['graph'].edges[e]['flow']


def get_excess_flow(path, mfd):

    if not path:
        return 0

    excess_flow = get_flow(path[0], mfd)
    for e in path[1:]:
        excess_flow = excess_flow + get_flow(e, mfd) - mfd['out_flow'][e[0]]

    return excess_flow


def find_right_maximal_extension_scan(mfd, paths, path, first, last):

    while last + 1 < len(path) and is_safe(mfd, len(paths), path, first, last + 2):
        last += 1

    return last


def find_right_maximal_extension_bin_search(mfd, paths, path, first, last):
    return bisect(range(last + 1, len(path)), 0, key=lambda i: 0 if is_safe(mfd, len(paths), path, first, i + 1) else 1) + last


def find_right_maximal_extension_exp_search(mfd, paths, path, first, last):

    exp = 0
    while last + 1 + exp < len(path) and is_safe(mfd, len(paths), path, first, last + 2 + exp):
        exp = 2 * exp + 1

    if exp < 3:
        return last + exp

    exp = int((exp-1)/2)

    return bisect(range(last + 2 + exp, min(last + 1 + 2*exp+1, len(path))), 0, key=lambda i: 0 if is_safe(mfd, len(paths), path, first, i + 1) else 1) + last + 1 + exp


def find_right_maximal_extension_repeated_exp_search(mfd, paths, path, first, last):
    return find_right_maximal_extension_repeated_exp_search_rec(mfd, paths, path, first, last, len(path))


def find_right_maximal_extension_repeated_exp_search_rec(mfd, paths, path, first, last, limit):

    exp = 0
    while last + 1 + exp < limit and is_safe(mfd, len(paths), path, first, last + 2 + exp):
        exp = 2 * exp + 1

    if exp < 3:
        return last + exp

    exp = int((exp - 1) / 2)

    return find_right_maximal_extension_repeated_exp_search_rec(mfd, paths, path, first, last + 1 + exp, min(last + 1 + 2 * exp + 1, len(path)))


def find_left_minimal_reduction_scan(mfd, paths, path, first, last, limit=None):

    first += 1
    while first < limit if limit else last + 1 and not is_safe(mfd, len(paths), path, first, last + 2):
        first += 1

    return first


def find_left_minimal_reduction_bin_search(mfd, paths, path, first, last, limit=None):

    return bisect(range(first + 1, limit if limit else last + 1), 0, key=lambda i: 1 if is_safe(mfd, len(paths), path, i, last + 2) else 0) + first + 1


def find_left_minimal_reduction_exp_search(mfd, paths, path, first, last, limit=None):

    exp = 0
    while first + 1 + exp < limit if limit else last + 1 and not is_safe(mfd, len(paths), path, first+1+exp, last + 2):
        exp = 2 * exp + 1

    if exp < 3:
        return first + 1 + exp

    exp = int((exp - 1) / 2)

    return bisect(range(first + 2 + exp, first + 1 + 2*exp+1), 0, key=lambda i: 1 if is_safe(mfd, len(paths), path, i, last + 2) else 0) + first + 2 + exp


def find_left_minimal_reduction_repeated_exp_search(mfd, paths, path, first, last, limit=None):
    return find_left_minimal_reduction_repeated_exp_search_rec(mfd, paths, path, first, last, limit if limit else last + 1)


def find_left_minimal_reduction_repeated_exp_search_rec(mfd, paths, path, first, last, limit):

    exp = 0
    while first + 1 + exp < limit and not is_safe(mfd, len(paths), path, first + 1 + exp, last + 2):
        exp = 2 * exp + 1

    if exp < 3:
        return first + 1 + exp

    exp = int((exp - 1) / 2)

    return find_left_minimal_reduction_repeated_exp_search_rec(mfd, paths, path, first + 1 + exp, last, first + 1 + 2 * exp + 1)


def compute_maximal_safe_paths(mfd):

    paths = mfd['solution']
    max_safe_paths = list()

    for path in paths:

        maximal_safe_paths = list()

        # two-finger algorithm
        # Invariant: path[first:last+1] (python notation) is safe
        first = 0
        last = 0

        while True:

            # Extending
            if len(path) - last > mfd['extension_strategy_threshold']:
                last = mfd['extension_strategy_large'](mfd, paths, path, first, last)
            else:
                last = mfd['extension_strategy_small'](mfd, paths, path, first, last)

            maximal_safe_paths.append((first, last + 1))
            if last == len(path) - 1:
                break

            # Reducing
            if last+1 - first > mfd['reduction_strategy_threshold']:
                first = mfd['reduction_strategy_large'](mfd, paths, path, first, last)
            else:
                first = mfd['reduction_strategy_small'](mfd, paths, path, first, last)

            last += 1

        max_safe_paths.append(maximal_safe_paths)

    return paths, max_safe_paths


def compute_maximal_safe_paths_using_excess_flow(mfd):

    paths = mfd['solution']
    max_safe_paths = list()

    for path in paths:

        maximal_safe_paths = list()

        # two-finger algorithm
        # Invariant: path[first:last+1] (python notation) is safe
        first = 0
        last = 0

        while True:

            excess_flow = get_excess_flow(path[first:last+1], mfd)
            # Extending
            while last+1 < len(path) and excess_flow - mfd['out_flow'][path[last][1]] + get_flow(path[last + 1], mfd) > 0:
                excess_flow = excess_flow - mfd['out_flow'][path[last][1]] + get_flow(path[last + 1], mfd)
                last += 1

            if len(path) - last > mfd['extension_strategy_threshold']:
                last = mfd['extension_strategy_large'](mfd, paths, path, first, last)
            else:
                last = mfd['extension_strategy_small'](mfd, paths, path, first, last)

            maximal_safe_paths.append((first, last + 1))

            if last == len(path) - 1:
                break

            # Reducing
            ef = get_flow(path[last+1], mfd)
            limit = last+1

            while limit > 0 and ef - mfd['out_flow'][path[limit][0]] + get_flow(path[limit - 1], mfd) > 0:
                ef = ef - mfd['out_flow'][path[limit][0]] + get_flow(path[limit - 1], mfd)
                limit -= 1

            if limit - first > mfd['reduction_strategy_threshold']:
                first = mfd['reduction_strategy_large'](mfd, paths, path, first, last, limit)
            else:
                first = mfd['reduction_strategy_small'](mfd, paths, path, first, last, limit)

            last += 1

        max_safe_paths.append(maximal_safe_paths)

    return paths, max_safe_paths


def get_out_tree(ngraph, v, processed, contracted):

    cont = {e: contracted[e] for e in contracted}

    # Get root of out_tree
    root = v
    while ngraph.in_degree(root) == 1:
        root = next(ngraph.predecessors(root))

    leaf_edges = list()
    edges = deque(ngraph.out_edges(root, keys=True))

    while edges:
        u, v, i = edges.popleft()
        processed[u] = True
        cont[u, v, i] = True
        if ngraph.in_degree(v) > 1 or ngraph.out_degree(v) == 0:
            leaf_edges.append((u, v, i))
        else:
            for e in ngraph.out_edges(v, keys=True):
                edges.append(e)

    # Filter out uncompressible edges
    for e in [(u, v, i) for (u, v, i) in leaf_edges if u == root]:
        cont[e] = contracted[e]
    leaf_edges = [(u, v, i) for (u, v, i) in leaf_edges if u != root]

    return root, leaf_edges, processed, cont


def get_in_tree(ngraph, v, processed, contracted):

    cont = {e: contracted[e] for e in contracted}

    # Get root of in_tree
    root = v
    while ngraph.out_degree(root) == 1:
        root = next(ngraph.successors(root))

    leaf_edges = list()
    edges = deque(ngraph.in_edges(root, keys=True))

    while edges:
        u, v, i = edges.popleft()
        processed[v] = True
        cont[u, v, i] = True
        if ngraph.out_degree(u) > 1 or ngraph.in_degree(u) == 0:
            leaf_edges.append((u, v, i))
        else:
            for e in ngraph.in_edges(u, keys=True):
                edges.append(e)

    # Filter out uncompressible edges
    for e in [(u, v, i) for (u, v, i) in leaf_edges if v == root]:
        cont[e] = contracted[e]
    leaf_edges = [(u, v, i) for (u, v, i) in leaf_edges if v != root]

    return root, leaf_edges, processed, cont


def get_out_trees(ngraph):

    out_trees = dict()
    contracted = {e: False for e in ngraph.edges(keys=True)}
    processed = {v: ngraph.in_degree(v) != 1 for v in ngraph.nodes}
    for v in ngraph.nodes:
        if not processed[v]:
            root, leaf_edges, processed, contracted = get_out_tree(ngraph, v, processed, contracted)
            if leaf_edges:
                out_trees[root] = leaf_edges
    return out_trees, contracted


def get_in_trees(ngraph):

    in_trees = dict()
    contracted = {e: False for e in ngraph.edges(keys=True)}
    processed = {v: ngraph.out_degree(v) != 1 for v in ngraph.nodes}
    for v in ngraph.nodes:
        if not processed[v]:
            root, leaf_edges, processed, contracted = get_in_tree(ngraph, v, processed, contracted)
            if leaf_edges:
                in_trees[root] = leaf_edges
    return in_trees, contracted


def get_y_to_v_contraction(ngraph):

    out_trees, contracted = get_out_trees(ngraph)

    mngraph_out_contraction = nx.MultiDiGraph()
    for u, v, i in ngraph.edges(keys=True):
        if not contracted[u, v, i]:
            mngraph_out_contraction.add_edge(u, v, flow=ngraph.edges[u, v, i]['flow'], pred=u)

    for root, leaf_edges in out_trees.items():
        for u, v, i in leaf_edges:
            mngraph_out_contraction.add_edge(root, v, flow=ngraph.edges[u, v, i]['flow'], pred=u)

    in_trees, contracted = get_in_trees(mngraph_out_contraction)

    mngraph_in_contraction = nx.MultiDiGraph()
    for u, v, i in mngraph_out_contraction.edges(keys=True):
        if not contracted[u, v, i]:
            mngraph_in_contraction.add_edge(u, v, flow=mngraph_out_contraction.edges[u, v, i]['flow'], succ=(v, i))

    for root, leaf_edges in in_trees.items():
        for u, v, i in leaf_edges:
            mngraph_in_contraction.add_edge(u, root, flow=mngraph_out_contraction.edges[u, v, i]['flow'], succ=(v, i))

    # Remove trivial_paths found as edges from source to sink in the contraction
    trivial_paths = list()
    edges = list(mngraph_in_contraction.edges(keys=True))
    for u, v, i in edges:
        if mngraph_in_contraction.in_degree(u) == 0 and mngraph_in_contraction.out_degree(v) == 0:
            trivial_paths.append(get_expanded_path([(u, v, i)], mngraph_in_contraction, ngraph, mngraph_out_contraction))
            mngraph_in_contraction.remove_edge(u, v, i)

    return mngraph_out_contraction, mngraph_in_contraction, trivial_paths


def get_expanded_path(path, graph, original_graph, out_contraction_graph):

    out_contraction_path = list()

    for u, v, i in path:
        root = v
        v, i = graph.edges[u, v, i]['succ']

        expanded_edge = list()
        while v != root:
            expanded_edge.append((u, v, i))
            u, v, i = list(out_contraction_graph.out_edges(v, keys=True))[0]

        expanded_edge.append((u, v, i))
        out_contraction_path += list(expanded_edge)

    original_path = list()

    for u, v, i in out_contraction_path:
        root = u
        u = out_contraction_graph.edges[u, v, i]['pred']
        expanded_edge = list()
        while u != root:
            expanded_edge.append((u, v, i))
            u, v, i = list(original_graph.in_edges(u, keys=True))[0]
        expanded_edge.append((u, v, i))
        original_path += list(reversed(expanded_edge))

    return original_path


def compute_graph_metadata(graph, use_excess_flow, use_y_to_v, strategy):

    # creation of NetworkX Graph
    ngraph = nx.MultiDiGraph()
    ngraph.add_weighted_edges_from(graph['edges'], weight='flow')
    if use_y_to_v:
        original = ngraph
        out_contraction, ngraph, trivial_paths = get_y_to_v_contraction(ngraph)
        mapping = (original, out_contraction)
        print(f'Y2V contraction reduced the number of edges from {len(original.edges)} to {len(ngraph.edges)}')

    # calculating source, sinks
    sources = [x for x in ngraph.nodes if ngraph.in_degree(x) == 0]
    sinks = [x for x in ngraph.nodes if ngraph.out_degree(x) == 0]

    if use_excess_flow:
        in_flow = {v: sum([f for u, v, f in ngraph.in_edges(v, data='flow')]) for v in ngraph.nodes}
        out_flow = {u: sum([f for u, v, f in ngraph.out_edges(u, data='flow')]) for u in ngraph.nodes}

    # definition of data
    return {
        'graph': ngraph,
        'sources': sources,
        'sinks': sinks,
        'max_flow_value': max(ngraph.edges(data='flow'), key=lambda e: e[-1])[-1] if len(ngraph.edges) > 0 else -1,
        'compute': compute_maximal_safe_paths_using_excess_flow if use_excess_flow else compute_maximal_safe_paths,
        'in_flow': in_flow if use_excess_flow else None,
        'out_flow': out_flow if use_excess_flow else None,
        'mapping': mapping if use_y_to_v else None,
        'trivial_paths': trivial_paths if use_y_to_v else None,
        'extension_strategy_small': find_right_maximal_extension_scan if 'extension_strategy_small' not in strategy else strategy['extension_strategy_small'],
        'extension_strategy_large': find_right_maximal_extension_scan if 'extension_strategy_large' not in strategy else strategy['extension_strategy_large'],
        'reduction_strategy_small': find_left_minimal_reduction_scan if 'reduction_strategy_small' not in strategy else strategy['reduction_strategy_small'],
        'reduction_strategy_large': find_left_minimal_reduction_scan if 'reduction_strategy_large' not in strategy else strategy['reduction_strategy_large'],
        'reduction_strategy_threshold': float('-inf') if 'reduction_strategy_threshold' not in strategy else strategy['reduction_strategy_threshold'],
        'extension_strategy_threshold': float('-inf') if 'extension_strategy_threshold' not in strategy else strategy['extension_strategy_threshold'],
    }


def solve_instances_safety(graphs, output_file, use_excess_flow, use_y_to_v, strategy):

    output = open(output_file, 'w+')
    output_counters = open(f'{output_file}.count', 'w+')

    for g, graph in enumerate(graphs):

        output.write(f'# graph {g}\n')
        output_counters.write(f'# graph {g}\n')

        if not graph['edges']:
            continue

        mfd = compute_graph_metadata(graph, use_excess_flow, use_y_to_v, strategy)
        global ilp_counter
        ilp_counter = 0

        if len(mfd['graph'].edges) > 0:

            mfd = mfd_algorithm(mfd)
            paths, max_safe_paths = mfd['compute'](mfd)

            for path, maximal_safe_paths in zip(paths, max_safe_paths):
                # print path
                for i, j in maximal_safe_paths:
                    max_path = path[i:j]
                    if use_y_to_v:
                        max_path = get_expanded_path(max_path, mfd['graph'], mfd['mapping'][0], mfd['mapping'][1])
                    output_maximal_safe_path(output, max_path)

        if use_y_to_v:
            for trivial_path in mfd['trivial_paths']:
                output_maximal_safe_path(output, trivial_path)

        output_counters.write(f'{ilp_counter}\n')

    output.close()
    output_counters.close()


def parse_strategy(strategy_str):

    if strategy_str == 'scan':
        return find_right_maximal_extension_scan, find_left_minimal_reduction_scan

    if strategy_str == 'bin_search':
        return find_right_maximal_extension_bin_search, find_left_minimal_reduction_bin_search

    if strategy_str == 'exp_search':
        return find_right_maximal_extension_exp_search, find_left_minimal_reduction_exp_search

    if strategy_str == 'rep_exp_search':
        return find_right_maximal_extension_repeated_exp_search, find_left_minimal_reduction_repeated_exp_search

    return None


def get_strategy(args):

    strategy = dict()

    if args.strategy_threshold:
        strategy['extension_strategy_threshold'] = args.strategy_threshold
        strategy['reduction_strategy_threshold'] = args.strategy_threshold

    if args.extension_strategy_threshold:
        strategy['extension_strategy_threshold'] = args.extension_strategy_threshold

    if args.reduction_strategy_threshold:
        strategy['reduction_strategy_threshold'] = args.reduction_strategy_threshold

    if parse_strategy(args.strategy):
        strategy['extension_strategy_small'] = parse_strategy(args.strategy)[0]
        strategy['extension_strategy_large'] = parse_strategy(args.strategy)[0]
        strategy['reduction_strategy_small'] = parse_strategy(args.strategy)[1]
        strategy['reduction_strategy_large'] = parse_strategy(args.strategy)[1]

    if parse_strategy(args.extension_strategy):
        strategy['extension_strategy_small'] = parse_strategy(args.extension_strategy)[0]
        strategy['extension_strategy_large'] = parse_strategy(args.extension_strategy)[0]

    if parse_strategy(args.reduction_strategy):
        strategy['reduction_strategy_small'] = parse_strategy(args.reduction_strategy)[1]
        strategy['reduction_strategy_large'] = parse_strategy(args.reduction_strategy)[1]

    if parse_strategy(args.extension_strategy_small):
        strategy['extension_strategy_small'] = parse_strategy(args.extension_strategy_small)[0]

    if parse_strategy(args.extension_strategy_large):
        strategy['extension_strategy_large'] = parse_strategy(args.extension_strategy_large)[0]

    if parse_strategy(args.reduction_strategy_small):
        strategy['reduction_strategy_small'] = parse_strategy(args.reduction_strategy_small)[1]

    if parse_strategy(args.reduction_strategy_large):
        strategy['reduction_strategy_large'] = parse_strategy(args.reduction_strategy_large)[1]

    return strategy


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

    parser.add_argument('-s', '--strategy', type=str, help='Strategy for extension and reduction of two-finger algorithm {scan, bin_search, exp_search, rep_exp_search}')
    parser.add_argument('-es', '--extension-strategy', type=str, help='Strategy for extension of two-finger algorithm {scan, bin_search, exp_search, rep_exp_search}')
    parser.add_argument('-rs', '--reduction-strategy', type=str, help='Strategy for reduction of two-finger algorithm {scan, bin_search, exp_search, rep_exp_search}')
    parser.add_argument('-ess', '--extension-strategy-small', type=str, help='Strategy for extension of two-finger algorithm {scan, bin_search, exp_search, rep_exp_search} when the search space is small (specified by threshold)')
    parser.add_argument('-esl', '--extension-strategy-large', type=str, help='Strategy for extension of two-finger algorithm {scan, bin_search, exp_search, rep_exp_search} when the search space is large (specified by threshold)')
    parser.add_argument('-rss', '--reduction-strategy-small', type=str, help='Strategy for reduction of two-finger algorithm {scan, bin_search, exp_search, rep_exp_search} when the search space is small (specified by threshold)')
    parser.add_argument('-rsl', '--reduction-strategy-large', type=str, help='Strategy for reduction of two-finger algorithm {scan, bin_search, exp_search, rep_exp_search} when the search space is large (specified by threshold)')
    parser.add_argument('-st', '--strategy-threshold', type=int, help='Search space threshold to switch from --<>-strategy-small to --<>-strategy-large')
    parser.add_argument('-est', '--extension-strategy-threshold', type=int, help='Search space threshold to switch from --extension-strategy-small to --extension-strategy-large')
    parser.add_argument('-rst', '--reduction-strategy-threshold', type=int, help='Search space threshold to switch from --reduction-strategy-small to --reduction-strategy-large')

    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-i', '--input', type=str, help='Input filename', required=True)
    requiredNamed.add_argument('-o', '--output', type=str, help='Output filename', required=True)

    args = parser.parse_args()

    threads = args.threads
    if threads == 0:
        threads = os.cpu_count()
    print(f'INFO: Using {threads} threads for the Gurobi solver')

    solve_instances_safety(read_input(args.input), args.output, args.use_excess_flow, args.use_y_to_v, get_strategy(args))
