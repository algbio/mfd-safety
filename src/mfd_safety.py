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
from copy import deepcopy


class TimeoutILP(Exception):
    pass


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

        # Look for subpaths in the graph
        try:
            subpaths_title_line_number = lines.index('subpaths')
        except ValueError:
            subpaths_title_line_number = len(lines)

        graph['n'] = int(lines[0])
        graph['edges'] = [get_edge(raw_e) for raw_e in lines[1:subpaths_title_line_number]]
        graph['subpaths'] = [[int(v) for v in raw_sub.split()[:-1]] for raw_sub in
                             lines[subpaths_title_line_number + 1:]]
        graph['subpaths'] = [[(sp[i], sp[i + 1], 0) for i in range(len(sp) - 1)] for sp in graph['subpaths']]

    finally:
        return graph


def read_input_graphs(graph_file):
    graphs_raw = open(graph_file, 'r').read().split('#')[1:]
    return [get_graph(raw_g) for raw_g in graphs_raw]


def read_input(graph_file):
    return read_input_graphs(graph_file)


def mfd_algorithm_repeated_exp_search(data, base_unfeasible, sequential_threshold):
    exp = 1
    while fd_fixed_size(data, base_unfeasible + exp)['message'] != 'solved':
        exp *= 2

    if exp <= 2 * sequential_threshold:
        data_copy = deepcopy(data)
        for i in range(base_unfeasible + int(exp / 2) + 1, base_unfeasible + exp):
            if fd_fixed_size(data, i)['message'] == 'solved':
                return i
        data.update(data_copy)
        return base_unfeasible + exp

    return mfd_algorithm_repeated_exp_search(data, base_unfeasible + int(exp / 2), sequential_threshold)


def mfd_algorithm(data):
    sequential_threshold = data['sequential_threshold']

    data['message'] = 'unsolved'
    if sequential_threshold == 0:
        for i in range(1, len(data['graph'].edges) + 1):
            if fd_fixed_size(data, i)['message'] == 'solved':
                return data
    else:
        mfd_algorithm_repeated_exp_search(data, 0, sequential_threshold)

    return data


def build_base_ilp_model(data, size):
    graph = data['graph']
    max_flow_value = data['max_flow_value']
    sources = data['sources']
    sinks = data['sinks']

    # create extra sets
    t = [(u, v, i, k) for (u, v, i) in graph.edges(keys=True) for k in range(size)]

    # Create a new model
    model = gp.Model('MFD')
    model.setParam('LogToConsole', 0)
    model.setParam('Threads', threads)

    # Create variables
    x = model.addVars(t, vtype=GRB.BINARY, name='x')
    w = model.addVars(range(size), vtype=GRB.INTEGER, name='w', lb=0)
    z = model.addVars(t, vtype=GRB.CONTINUOUS, name='z', lb=0)

    # flow conservation
    for k in range(size):
        for v in graph.nodes:
            if v in sources:
                model.addConstr(sum(x[v, w, i, k] for _, w, i in graph.out_edges(v, keys=True)) == 1)
            if v in sinks:
                model.addConstr(sum(x[u, v, i, k] for u, _, i in graph.in_edges(v, keys=True)) == 1)
            if v not in sources and v not in sinks:
                model.addConstr(sum(x[v, w, i, k] for _, w, i in graph.out_edges(v, keys=True)) - sum(
                    x[u, v, i, k] for u, _, i in graph.in_edges(v, keys=True)) == 0)

    # flow balance
    for (u, v, i, f) in graph.edges(keys=True, data='flow'):
        model.addConstr(f == sum(z[u, v, i, k] for k in range(size)))

    # linearization
    for (u, v, i) in graph.edges(keys=True):
        for k in range(size):
            model.addConstr(z[u, v, i, k] <= max_flow_value * x[u, v, i, k])
            model.addConstr(w[k] - (1 - x[u, v, i, k]) * max_flow_value <= z[u, v, i, k])
            model.addConstr(z[u, v, i, k] <= w[k])

    subpaths = data['subpaths']
    if subpaths:
        t = [(s, k) for s in range(len(subpaths)) for k in range(size)]
        rho = model.addVars(t, vtype=GRB.BINARY, name="rho")

        # subpath constraints
        for s, path in enumerate(subpaths):
            for k in range(size):
                model.addConstr(sum(x[u, v, i, k] for (u, v, i) in path) >= len(path) * rho[s, k])
                # guarantee that all edges of a subpath are selected

            model.addConstr(sum(rho[s, k] for k in range(size)) >= 1)
            # subpath is present in at least one of the size paths

    return model, x, w, z


def build_ilp_model_avoiding_multiple_paths(data, size, paths):
    model, x, _, _ = build_base_ilp_model(data, size)

    r = model.addVars(range(len(paths)), vtype=GRB.BINARY, name="r")
    model.setObjective(sum(r[p] for p in range(len(paths))), GRB.MAXIMIZE)

    # group testing constraint
    for p, path in enumerate(paths):
        for k in range(size):
            model.addConstr(sum(x[u, v, i, k] for (u, v, i) in path) <= len(path) - r[p])

    return model


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
            w_sol[k] = round(model.getVarByName(f'w[{k}]').x)
        for (u, v, i, k) in T:
            if round(model.getVarByName(f'x[{u},{v},{i},{k}]').x) == 1:
                paths[k].append((u, v, i))
        for k in range(len(paths)):
            paths[k] = sorted(paths[k])

        data['weights'], data['solution'] = w_sol, paths

    return data


def update_status(data, model):
    global ilp_counter
    ilp_counter += 1

    if model.status == GRB.OPTIMAL:
        data['message'] = 'solved'
        data['runtime'] = model.Runtime

    if model.status == GRB.INFEASIBLE:
        data['message'] = 'unsolved'
        data['runtime'] = 0

    if ilp_time_budget:
        global time_budget
        time_budget -= model.Runtime

        if model.status == GRB.TIME_LIMIT:
            raise TimeoutILP()

    return data


def fd_fixed_size_forbidding_path(data, size, path):
    if ilp_time_budget and time_budget < 0:
        raise TimeoutILP()

    # calculate a flow decomposition into size paths avoiding path
    try:
        # Create a new model
        model = build_ilp_model_avoiding_path(data, size, path)
        if ilp_time_budget:
            model.setParam('TimeLimit', time_budget)
        model.optimize()

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
    if ilp_time_budget and time_budget < 0:
        raise TimeoutILP()

    # calculate a flow decomposition into size paths
    try:
        # Create a new model
        model, _, _, _ = build_base_ilp_model(data, size)

        if ilp_time_budget:
            model.setParam('TimeLimit', time_budget)

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
    return bisect(range(last + 1, len(path)), 0,
                  key=lambda i: 0 if is_safe(mfd, len(paths), path, first, i + 1) else 1) + last


def find_right_maximal_extension_exp_search(mfd, paths, path, first, last):
    exp = 0
    while last + 1 + exp < len(path) and is_safe(mfd, len(paths), path, first, last + 2 + exp):
        exp = 2 * exp + 1

    if exp < 3:
        return last + exp

    exp = int((exp - 1) / 2)

    return bisect(range(last + 2 + exp, min(last + 1 + 2 * exp + 1, len(path))), 0,
                  key=lambda i: 0 if is_safe(mfd, len(paths), path, first, i + 1) else 1) + last + 1 + exp


def find_right_maximal_extension_repeated_exp_search(mfd, paths, path, first, last):
    return find_right_maximal_extension_repeated_exp_search_rec(mfd, paths, path, first, last, len(path))


def find_right_maximal_extension_repeated_exp_search_rec(mfd, paths, path, first, last, limit):
    exp = 0
    while last + 1 + exp < limit and is_safe(mfd, len(paths), path, first, last + 2 + exp):
        exp = 2 * exp + 1

    if exp < 3:
        return last + exp

    exp = int((exp - 1) / 2)

    return find_right_maximal_extension_repeated_exp_search_rec(mfd, paths, path, first, last + 1 + exp,
                                                                min(last + 1 + 2 * exp + 1, len(path)))


def find_left_minimal_reduction_scan(mfd, paths, path, first, last, limit):
    first += 1
    while first < limit and not is_safe(mfd, len(paths), path, first, last + 2):
        first += 1

    return first


def find_left_minimal_reduction_bin_search(mfd, paths, path, first, last, limit):
    return bisect(range(first + 1, limit), 0,
                  key=lambda i: 1 if is_safe(mfd, len(paths), path, i, last + 2) else 0) + first + 1


def find_left_minimal_reduction_exp_search(mfd, paths, path, first, last, limit):
    exp = 0
    while first + 1 + exp < limit and not is_safe(mfd, len(paths), path, first + 1 + exp, last + 2):
        exp = 2 * exp + 1

    if exp < 3:
        return first + 1 + exp

    exp = int((exp - 1) / 2)

    return bisect(range(first + 2 + exp, first + 1 + 2 * exp + 1), 0,
                  key=lambda i: 1 if is_safe(mfd, len(paths), path, i, last + 2) else 0) + first + 2 + exp


def find_left_minimal_reduction_repeated_exp_search(mfd, paths, path, first, last, limit):
    return find_left_minimal_reduction_repeated_exp_search_rec(mfd, paths, path, first, last, limit)


def find_left_minimal_reduction_repeated_exp_search_rec(mfd, paths, path, first, last, limit):
    exp = 0
    while first + 1 + exp < limit and not is_safe(mfd, len(paths), path, first + 1 + exp, last + 2):
        exp = 2 * exp + 1

    if exp < 3:
        return first + 1 + exp

    exp = int((exp - 1) / 2)

    return find_left_minimal_reduction_repeated_exp_search_rec(mfd, paths, path, first + 1 + exp, last,
                                                               first + 1 + 2 * exp + 1)


def get_unsafe_unknown(model, test):
    if model.status == GRB.OPTIMAL:
        return [test for p, test in enumerate(test) if round(model.getVarByName(f'r[{p}]').x) == 1], [test for p, test
                                                                                                      in enumerate(test)
                                                                                                      if round(
                model.getVarByName(f'r[{p}]').x) == 0]

    return list(), test


def get_trivially_safe_rest(mfd, k, paths, test, use_excess_flow):
    if use_excess_flow:
        fd_safe = [get_excess_flow(paths[p][i:j], mfd) > 0 for p, i, j in test]
        return [(p, i, j) for p, i, j in test if fd_safe[p]], [(p, i, j) for p, i, j in test if not fd_safe[p]]
    return [(p, i, j) for p, i, j in test if i == j + 1], [(p, i, j) for p, i, j in test if i != j + 1]


def get_safe_unsafe(mfd, k, paths, test):
    """Takes in a set of paths and returns those paths that are safe and those that are unsafe.
    paths is a list of paths of an mfd.
    test is a list of indices (lists of size 3), where each element p,i,j represents the path paths[p][i:j].
    Implements Algorithm 3 in the paper."""

    if ilp_time_budget and time_budget < 0:
        raise TimeoutILP()

    model = build_ilp_model_avoiding_multiple_paths(mfd, k, [paths[p][i:j] for p, i, j in test])

    if ilp_time_budget:
        model.setParam('TimeLimit', time_budget)
    model.optimize()

    mfd = update_status(mfd, model)

    unsafe, unknown = get_unsafe_unknown(model, test)

    if not unsafe:
        return unknown, unsafe
    else:
        s, u = get_safe_unsafe(mfd, k, paths, unknown)
        return s, unsafe + u


def get_extensions(safe_paths, paths):
    return list(set([(p, i - 1, j) for p, i, j in safe_paths if i >= 0] + [(p, i, j + 1) for p, i, j in safe_paths if
                                                                           j < len(paths[p])]))


def all_maximal_safe_paths_bottom_up(mfd, safe_paths, max_safe, use_excess_flow=False):
    """
    Use group testing in a bottom up fashion.
    Implements algorithm x from paper.
    """
    if not safe_paths:
        return

    paths = mfd['solution']
    test = get_extensions(safe_paths, paths)

    t_safe, test = get_trivially_safe_rest(mfd, len(paths), paths, test, use_excess_flow)
    safe, unsafe = get_safe_unsafe(mfd, len(paths), paths, test)
    safe += t_safe

    for p, i, j in safe_paths:
        right_maximal = j == len(paths[p]) or (p, i, j + 1) in unsafe
        left_maximal = i == 0 or (p, i - 1, j) in unsafe
        if right_maximal and left_maximal:
            max_safe.append((p, i, j))

    all_maximal_safe_paths_bottom_up(mfd, safe, max_safe)


def get_reductions(unsafe_paths, paths):
    return list(
        set([(p, i + 1, j) for p, i, j in unsafe_paths if j == len(paths[p]) or (p, i + 1, j + 1) in unsafe_paths] + [
            (p, i, j - 1) for p, i, j in unsafe_paths if (i >= 0) and (i == 0 or (i - 1, j - 1, p) in unsafe_paths)]))


def all_maximal_safe_paths_top_down(mfd, unsafe_paths, max_safe, use_excess_flow=False):
    """
    Use group testing in a top down fashion.
    Implements algorithm 6 from paper.
    """
    if not unsafe_paths:
        return

    paths = mfd['solution']
    test = get_reductions(unsafe_paths, paths)

    t_safe, test = get_trivially_safe_rest(mfd, len(paths), paths, test, use_excess_flow)
    safe, unsafe = get_safe_unsafe(mfd, len(paths), paths, test)
    safe += t_safe

    for p, i, j in safe:
        max_safe.append((p, i, j))

    all_maximal_safe_paths_top_down(mfd, unsafe, max_safe)


def compute_maximal_safe_paths_using_group_bottom_up(mfd):
    paths = mfd['solution']

    max_safe = list()
    all_maximal_safe_paths_bottom_up(mfd, [(p, i, i + 1) for p, path in enumerate(paths) for i in range(len(path))],
                                     max_safe, use_excess_flow=mfd['use_excess_flow'])

    maximal_safe_paths = [[] for path in paths]
    for p, i, j in max_safe:
        maximal_safe_paths[p].append((i, j))

    return paths, maximal_safe_paths


def compute_maximal_safe_paths_using_group_top_down(mfd):
    paths = mfd['solution']

    max_safe = list()
    all_maximal_safe_paths_top_down(mfd, [(p, -1, len(path)) for p, path in enumerate(paths)], max_safe,
                                    use_excess_flow=mfd['use_excess_flow'])

    maximal_safe_paths = [[] for path in paths]
    for p, i, j in max_safe:
        maximal_safe_paths[p].append((i, j))

    return paths, maximal_safe_paths


def compute_maximal_safe_paths_using_two_finger(mfd):
    paths = mfd['solution']
    use_excess_flow = mfd['use_excess_flow']

    max_safe_paths = list()

    for path in paths:

        maximal_safe_paths = list()

        first = 0
        last = 0

        while True:

            # Extending
            if use_excess_flow:
                excess_flow = get_excess_flow(path[first:last + 1], mfd)

                while last + 1 < len(path) and excess_flow - mfd['out_flow'][path[last][1]] + get_flow(path[last + 1],
                                                                                                       mfd) > 0:
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
            if use_excess_flow:
                ef = get_flow(path[last + 1], mfd)
                limit = last + 1

                while limit > 0 and ef - mfd['out_flow'][path[limit][0]] + get_flow(path[limit - 1], mfd) > 0:
                    ef = ef - mfd['out_flow'][path[limit][0]] + get_flow(path[limit - 1], mfd)
                    limit -= 1
            else:
                limit = last + 1

            if limit - first > mfd['reduction_strategy_threshold']:
                first = mfd['reduction_strategy_large'](mfd, paths, path, first, last, limit)
            else:
                first = mfd['reduction_strategy_small'](mfd, paths, path, first, last, limit)

            last += 1

        max_safe_paths.append(maximal_safe_paths)

    return paths, max_safe_paths


def extract_weighted_path(sources, graph):
    path = list()
    for s in sources:
        first_edges = [e for e in graph.out_edges(s, keys=True) if graph.edges[e]['remaining_flow'] > 0]
        if first_edges:
            path.append(first_edges[0])
            while [e for e in graph.out_edges(path[-1][1], keys=True) if graph.edges[e]['remaining_flow'] > 0]:
                path.append(
                    [e for e in graph.out_edges(path[-1][1], keys=True) if graph.edges[e]['remaining_flow'] > 0][0])
            break

    return path, min([graph.edges[e]['remaining_flow'] for e in path])


def decompose_flow(mfd):
    weights, paths = list(), list()

    sources, sinks = mfd['sources'], mfd['sinks']
    graph = mfd['graph']

    remaining_flow = sum([mfd['out_flow'][s] for s in sources])

    edges = list(graph.edges(keys=True))
    for e in edges:
        graph.edges[e]['remaining_flow'] = graph.edges[e]['flow']

    while remaining_flow > 0:

        path, weight = extract_weighted_path(sources, graph)
        weights.append(weight)
        paths.append(path)

        for e in path:
            graph.edges[e]['remaining_flow'] -= weight
        remaining_flow -= weight

    mfd['solution'], mfd['weights'] = paths, weights

    return mfd


def compute_maximal_safe_paths_for_flow_decompositions(mfd):
    if 'in_flow' not in mfd:
        mfd['in_flow'] = {v: sum([f for u, v, f in mfd['graph'].in_edges(v, data='flow')]) for v in mfd['graph'].nodes}
        mfd['out_flow'] = {u: sum([f for u, v, f in mfd['graph'].out_edges(u, data='flow')]) for u in
                           mfd['graph'].nodes}

    if not mfd['solution']:
        mfd = decompose_flow(mfd)

    paths = mfd['solution']
    max_safe_paths = list()

    for path in paths:

        maximal_safe_paths = list()

        # two-finger algorithm
        # Invariant: path[first:last+1] (python notation) is safe
        first = 0
        last = 0
        excess_flow = get_excess_flow(path[first:last + 1], mfd)

        while True:

            # Extending
            while last + 1 < len(path) and excess_flow - mfd['out_flow'][path[last][1]] + get_flow(path[last + 1],
                                                                                                   mfd) > 0:
                excess_flow = excess_flow - mfd['out_flow'][path[last][1]] + get_flow(path[last + 1], mfd)
                last += 1

            maximal_safe_paths.append((first, last + 1))

            if last == len(path) - 1:
                break

            excess_flow = excess_flow - mfd['out_flow'][path[last][1]] + get_flow(path[last + 1], mfd)
            last += 1

            # Reducing
            while first < last and excess_flow < 0:
                excess_flow = excess_flow + mfd['in_flow'][path[first][1]] + get_flow(path[first], mfd)
                first += 1

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
            trivial_paths.append(
                get_expanded_path([(u, v, i)], mngraph_in_contraction, ngraph, mngraph_out_contraction))
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


def get_contracted_path(path, graph, original_graph, out_contraction_graph):
    # This code assumes that path starts and ends at a vertex in graph
    out_contraction_path = list()
    root = None
    for e in path:
        if e in out_contraction_graph.edges(keys=True):
            out_contraction_path.append(e)
        else:
            if root is None:
                root = e[0]

            if (root, e[1], e[0]) in out_contraction_graph.edges(data='pred'):
                candidates = [
                    ed[:-1] for ed in out_contraction_graph.edges(keys=True, data='pred')
                    if ed[0] == root and ed[1] == e[1] and ed[3] == e[0]
                ]

                # Here it assumes that the last one added is the corresponding contracted edge
                out_contraction_path.append(max(candidates, key=lambda ed: ed[2]))

                root = None

    reversed_contracted_path = list()
    root = None
    for e in out_contraction_path[::-1]:
        if e in graph.edges(keys=True):
            reversed_contracted_path.append(e)
        else:
            if root is None:
                root = e[1]

            if (e[0], root, (e[1], e[2])) in graph.edges(data='succ'):
                candidates = [
                    ed[:-1] for ed in graph.edges(keys=True, data='succ')
                    if ed[0] == e[0] and ed[1] == root and ed[3] == (e[1], e[2])
                ]

                # Here it assumes that the last one added is the corresponding contracted edge
                reversed_contracted_path.append(max(candidates, key=lambda ed: ed[2]))

                root = None

    return reversed_contracted_path[::-1]


def compute_graph_metadata(graph, use_excess_flow=False, use_y_to_v=False,
                           use_group_top_down=False, use_group_bottom_up=False,
                           sequential_threshold=0, strategy=dict()):
    # creation of NetworkX Graph
    ngraph = nx.MultiDiGraph()
    ngraph.add_weighted_edges_from(graph['edges'], weight='flow')
    subpaths = graph['subpaths']
    if use_y_to_v:
        original = ngraph
        out_contraction, ngraph, trivial_paths = get_y_to_v_contraction(ngraph)
        mapping = (original, out_contraction)
        subpaths = [get_contracted_path(path, ngraph, original, out_contraction) for path in subpaths]

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
        'compute': compute_maximal_safe_paths_using_group_bottom_up if use_group_bottom_up else compute_maximal_safe_paths_using_group_top_down if use_group_top_down else compute_maximal_safe_paths_using_two_finger,
        'in_flow': in_flow if use_excess_flow else None,
        'out_flow': out_flow if use_excess_flow else None,
        'mapping': mapping if use_y_to_v else None,
        'trivial_paths': trivial_paths if use_y_to_v else list(),
        'extension_strategy_small': find_right_maximal_extension_scan if 'extension_strategy_small' not in strategy else
        strategy['extension_strategy_small'],
        'extension_strategy_large': find_right_maximal_extension_scan if 'extension_strategy_large' not in strategy else
        strategy['extension_strategy_large'],
        'reduction_strategy_small': find_left_minimal_reduction_scan if 'reduction_strategy_small' not in strategy else
        strategy['reduction_strategy_small'],
        'reduction_strategy_large': find_left_minimal_reduction_scan if 'reduction_strategy_large' not in strategy else
        strategy['reduction_strategy_large'],
        'reduction_strategy_threshold': float('-inf') if 'reduction_strategy_threshold' not in strategy else strategy[
            'reduction_strategy_threshold'],
        'extension_strategy_threshold': float('-inf') if 'extension_strategy_threshold' not in strategy else strategy[
            'extension_strategy_threshold'],
        'use_excess_flow': use_excess_flow,
        'sequential_threshold': sequential_threshold,
        'subpaths': subpaths,
    }


def solve_instances_safety(graphs, output_file, output_stats=False, use_excess_flow=False,
                           use_y_to_v=False, use_group_top_down=False, use_group_bottom_up=False,
                           sequential_threshold=0, strategy=dict()):
    output = open(output_file, 'w+')
    if output_stats:
        stats = open(f'{output_file}.stats', 'w+')

    for g, graph in enumerate(graphs):

        output.write(f'# graph {g}\n')
        if output_stats:
            stats.write(f'# graph {g}\n')

        if not graph['edges']:
            continue

        mfd = compute_graph_metadata(graph, use_excess_flow, use_y_to_v,
                                     use_group_top_down, use_group_bottom_up,
                                     sequential_threshold, strategy)

        if output_stats and use_y_to_v:
            stats.write(f'Y2V: {len(mfd["graph"].edges)}/{len(mfd["mapping"][0].edges)}\n')

        global time_budget
        time_budget = ilp_time_budget

        global ilp_counter
        ilp_counter = 0

        if output_stats:
            is_unique_decomposition = True

        if len(mfd['graph'].edges) > 0:

            try:
                mfd = mfd_algorithm(mfd)
                paths, max_safe_paths = mfd['compute'](mfd)
                if output_stats:
                    stats.write('timeout: 0\n')

            except TimeoutILP:

                paths, max_safe_paths = compute_maximal_safe_paths_for_flow_decompositions(mfd)
                if output_stats:
                    stats.write('timeout: 1\n')

            for path, maximal_safe_paths in zip(paths, max_safe_paths):
                # print path
                for i, j in maximal_safe_paths:
                    max_path = path[i:j]
                    if use_y_to_v:
                        max_path = get_expanded_path(max_path, mfd['graph'], mfd['mapping'][0], mfd['mapping'][1])
                    output_maximal_safe_path(output, max_path)

                if output_stats and is_unique_decomposition:
                    is_unique_decomposition = len(maximal_safe_paths) == 1
                    if is_unique_decomposition:
                        max_path = maximal_safe_paths[0]
                        is_unique_decomposition = max_path == (0, len(path))

        else:
            if output_stats:
                stats.write('timeout: 0\n')

        if use_y_to_v:
            for trivial_path in mfd['trivial_paths']:
                output_maximal_safe_path(output, trivial_path)

        if output_stats:
            stats.write(
                f'Decomposition size: {(len(mfd["trivial_paths"]) if "trivial_paths" in mfd else 0) + (len(mfd["solution"]) if "solution" in mfd else 0)}\n')
            stats.write(f'ILP count: {ilp_counter}\n')
            if ilp_time_budget:
                stats.write(f'ILP time used: {ilp_time_budget - time_budget}/{ilp_time_budget}\n')

            if is_unique_decomposition:
                stats.write(f'Unique decomposition: 1\n')
            else:
                stats.write(f'Unique decomposition: 0\n')

    output.close()
    if output_stats:
        stats.close()


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
    parser.add_argument('-stats', '--output-stats', action='store_true', help='Output stats to file <output>.stats')
    parser.add_argument('-wt', '--weighttype', type=str, default='int+',
                        help='Type of path weights (default int+):\n   int+ (positive non-zero ints), \n   float+ (positive non-zero floats).')
    parser.add_argument('-t', '--threads', type=int, default=0,
                        help='Number of threads to use for the Gurobi solver; use 0 for all threads (default 0).')
    parser.add_argument('-seqt', '--sequential-threshold', type=int, default=0,
                        help='A repeated exponential search is performed to find the minimum flow decomposition, this parameter specifies the universe size at which a sequencial search is performed instead; use 0 to only perform sequential search (default 0).')

    parser.add_argument('-ilptb', '--ilp-time-budget', type=float,
                        help='Maximum time (in seconds) that the ilp solver is allowed to take when computing safe paths for one graph')

    parser.add_argument('-uef', '--use-excess-flow', action='store_true',
                        help='Use excess flow of a path to save ILP calls')
    parser.add_argument('-ugtd', '--use-group-top-down', action='store_true',
                        help='Use top down group testing')
    parser.add_argument('-ugbu', '--use-group-bottom-up', action='store_true',
                        help='Use bottom up group testing')
    parser.add_argument('-uy2v', '--use-y-to-v', action='store_true', help='Use Y to V contraction of the input graphs')

    parser.add_argument('-s', '--strategy', type=str,
                        help='Strategy for extension and reduction of two-finger algorithm {scan, bin_search, exp_search, rep_exp_search}')
    parser.add_argument('-es', '--extension-strategy', type=str,
                        help='Strategy for extension of two-finger algorithm {scan, bin_search, exp_search, rep_exp_search}')
    parser.add_argument('-rs', '--reduction-strategy', type=str,
                        help='Strategy for reduction of two-finger algorithm {scan, bin_search, exp_search, rep_exp_search}')
    parser.add_argument('-ess', '--extension-strategy-small', type=str,
                        help='Strategy for extension of two-finger algorithm {scan, bin_search, exp_search, rep_exp_search} when the search space is small (specified by threshold)')
    parser.add_argument('-esl', '--extension-strategy-large', type=str,
                        help='Strategy for extension of two-finger algorithm {scan, bin_search, exp_search, rep_exp_search} when the search space is large (specified by threshold)')
    parser.add_argument('-rss', '--reduction-strategy-small', type=str,
                        help='Strategy for reduction of two-finger algorithm {scan, bin_search, exp_search, rep_exp_search} when the search space is small (specified by threshold)')
    parser.add_argument('-rsl', '--reduction-strategy-large', type=str,
                        help='Strategy for reduction of two-finger algorithm {scan, bin_search, exp_search, rep_exp_search} when the search space is large (specified by threshold)')
    parser.add_argument('-st', '--strategy-threshold', type=int,
                        help='Search space threshold to switch from --<>-strategy-small to --<>-strategy-large')
    parser.add_argument('-est', '--extension-strategy-threshold', type=int,
                        help='Search space threshold to switch from --extension-strategy-small to --extension-strategy-large')
    parser.add_argument('-rst', '--reduction-strategy-threshold', type=int,
                        help='Search space threshold to switch from --reduction-strategy-small to --reduction-strategy-large')

    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-i', '--input', type=str, help='Input filename', required=True)
    requiredNamed.add_argument('-o', '--output', type=str, help='Output filename', required=True)

    args = parser.parse_args()

    threads = args.threads
    if threads == 0:
        threads = os.cpu_count()
    print(f'INFO: Using {threads} threads for the Gurobi solver')

    ilp_counter = 0
    ilp_time_budget = args.ilp_time_budget
    time_budget = args.ilp_time_budget
    solve_instances_safety(read_input(args.input), args.output, args.output_stats,
                           args.use_excess_flow, args.use_y_to_v,
                           args.use_group_top_down, args.use_group_bottom_up,
                           args.sequential_threshold, get_strategy(args))
