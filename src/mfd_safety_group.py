#!/usr/bin/env python
# coding: utf-8

import os
import sys
import argparse
import gurobipy as gp
from gurobipy import GRB
from mfd_safety import read_input_graphs, solve_instances, build_base_ilp_model


def get_path(raw_path):

    parts = raw_path.split()
    return [(int(parts[i]), int(parts[i + 1])) for i in range(len(parts) - 1)]


def get_group_test(raw_group_test):

    group_test = {
        'n': 0,
        'paths': list()
    }

    try:
        lines = raw_group_test.split('\n')[1:]
        if not lines[-1]:
            lines = lines[:-1]
        group_test['n'], group_test['paths'] = int(lines[0]), [get_path(raw_p) for raw_p in lines[1:]]

    finally:
        return group_test


def read_group_tests(group_testing_file):

    group_tests_raw = open(group_testing_file, 'r').read().split('#')[1:]
    return [get_group_test(raw_gt) for raw_gt in group_tests_raw]


def read_input(graph_file, group_testing_file):

    return read_input_graphs(graph_file), read_group_tests(group_testing_file)


def build_ilp_model_avoiding_multiple_paths(data, size, paths):

    model, x, _, _ = build_base_ilp_model(data, size)

    R = range(len(paths))
    r = model.addVars(R, vtype=GRB.BINARY, name="r")
    model.setObjective(sum(r[p] for p in range(len(paths))), GRB.MAXIMIZE)

    # group testing constraint
    for p, path in enumerate(paths):
        for k in range(size):
            model.addConstr(sum(x[i, j, k] for (i, j) in path) <= len(path) - r[p])

    return model


def get_unsafe_paths(model, paths):

    if model.status == GRB.OPTIMAL:
        return [path for p, path in enumerate(paths) if model.getVarByName(f'r[{p}]').x == 1]

    return list()


# renamed this function
def solve_group_testing_example(graphs, group_tests, output_file):

    mfds = solve_instances(graphs)

    output = open(output_file, 'w+')

    for g, (mfd, group_test) in enumerate(zip(mfds, group_tests)):

        output.write(f"# graph {g}\n")
        print("Paths we are group testing:")
        print(group_test["paths"])

        try:
            model = build_ilp_model_avoiding_multiple_paths(mfd, len(mfd['solution']), group_test['paths'])
            model.optimize()

            for unsafe_path in get_unsafe_paths(model, group_test['paths']):
                print(f"unsafe: {unsafe_path}")
                output.write(' '.join(map(str, sorted(set([item for t in unsafe_path for item in t])))))
                output.write('\n')

        except gp.GurobiError as e:
            print(f'Error code {e.errno}: {e}', file=sys.stderr)

        except AttributeError:
            print('Encountered an attribute error', file=sys.stderr)

    output.close()


def group_test_paths(paths, mfd):
    """
    Input: a set of paths in a graph, and the graph data (mfd variable)
    Output: the subset of the paths that are safe for MFD.
    This implements algorithm 3 from the MFD safety with ILP paper.
    """

    model = build_ilp_model_avoiding_multiple_paths(mfd,
                                                    len(mfd['solution']),
                                                    paths)
    model.optimize()

    N = get_unsafe_paths(model, paths)
    if not N:
        return paths
    else:
        # call on paths \ N
        return group_test_paths(list(set(paths).difference(set(N))), mfd)


def left_extensions(core):
    to_return = []
    for i, path in enumerate(core["paths"]):
        for (left, right) in core[i]:
            if left != 0:
                to_return.append(tuple(path[left - 1:right + 1]))  # right not inclusive
    return to_return


def right_extensions(core):
    to_return = []
    for i, path in enumerate(core["paths"]):
        for (left, right) in core[i]:
            if right != len(path) - 1:
                to_return.append(tuple(path[left:right + 2]))  # right not inclusive
    return to_return


def update_core(safe_paths, core):
    """
    Move the left and right indices for all core subpaths based on what was
    safe.
    """
    for i, path in enumerate(core["paths"]):
        new_core = []
        for (left, right) in core[i]:
            if left != 0:
                # check left extension
                subpath = tuple(path[left - 1:right + 1])
                if subpath in safe_paths:
                    new_core.append((left-1, right))
            if right != len(path) - 1:
                # check right extension
                subpath = tuple(path[left: right + 2])
                if subpath in safe_paths:
                    new_core.append((left, right + 1))
        core[i] = sorted(list(set(new_core)))


def get_maximal_safe(core, mfd, max_safe_paths):
    """
    Input: a core set of paths and graph data (mfd). If the core set is a safe
    core for the graph, then this returns all maximal safe paths.
    Output: maximal safe paths generated from the core set.
    This implements algorithm 4 from the MFD safety with ILP paper.
    """
    paths = list(set(left_extensions(core) + right_extensions(core)))
    safe_paths = group_test_paths(paths, mfd)
    # if the left extension is not safe and the right extension is not safe,
    # then the original path was safe
    for i, path in enumerate(core["paths"]):
        # print(f"checking for path {path}")
        for (left, right) in core[i]:
            # print(f" checking ({left}, {right})")
            # print(f"  left extension is {path[left - 1: right + 1]}")
            # print(f"  right extension is {path[left: right + 2]}")
            left_unsafe = tuple(path[left - 1: right + 1]) not in safe_paths
            right_unsafe = tuple(path[left: right + 2]) not in safe_paths
            if left_unsafe and right_unsafe:
                max_safe_paths.append(path[left:right + 1])
    update_core(safe_paths, core)
    # print("new core is:")
    # for i, p in enumerate(core["paths"]):
    #     print(core[i])

    # count how many core sets are left
    total_core = 0
    for i, path in enumerate(core["paths"]):
        total_core += len(core[i])

    if total_core > 0:
        get_maximal_safe(core, mfd, max_safe_paths)


def solve_group_testing(graphs, output_file):
    # goal for now: do a group test for each path with all length 2 subpaths.

    mfds = solve_instances(graphs)

    for g, mfd in enumerate(mfds):

        # for each path in mfd, test all len 2 subpaths
        paths = mfd['solution']
        # core is a dictionary associating a bunch of subpaths (identified by
        # their l and r) to each path. to start, core is just length 1
        # subpaths (single edges).
        # e.g., if the 0th path is [(0,1),(1,11),(11,3),(3,12)], then
        # core[0] = [(0,0), (1,1), (2,2), (3,3)]
        core = dict()
        core["paths"] = paths
        for i, path in enumerate(paths):
            # core[i] = all length 1 path indices for ith path
            core[i] = list(zip(list(range(len(path))),
                               list(range(len(path)))))

        max_safe_paths = []
        get_maximal_safe(core, mfd, max_safe_paths)
        for x in max_safe_paths:
            print(x)
        """
        for path in paths:
            len2paths = list(zip(path, path[1:]))
            try:
                safe_paths = group_test_paths(len2paths, mfd)
            except gp.GurobiError as e:
                print(f'Error code {e.errno}: {e}', file=sys.stderr)

            except AttributeError:
                print('Encountered an attribute error', file=sys.stderr)
                """


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="""
        Runs mfd-safety group testing for an input set of paths.
        This script uses the Gurobi ILP solver.
        """,
        formatter_class=argparse.RawTextHelpFormatter
        )
    parser.add_argument('-wt', '--weighttype', type=str, default='int+', help='Type of path weights (default int+):\n   int+ (positive non-zero ints), \n   float+ (positive non-zero floats).')
    parser.add_argument('-t', '--threads', type=int, default=0, help='Number of threads to use for the Gurobi solver; use 0 for all threads (default 0).')

    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-i', '--input', type=str, help='Input filename', required=True)
    requiredNamed.add_argument('-g', '--group-test', type=str, help='Input (group tests) filename', required=False)
    requiredNamed.add_argument('-o', '--output', type=str, help='Output filename', required=True)

    args = parser.parse_args()

    threads = args.threads
    if threads == 0:
        threads = os.cpu_count()
    print(f"INFO: Using {threads} threads for the Gurobi solver")

    # solve_group_testing(*read_input(args.input, args.group_test), args.output)
    solve_group_testing(read_input_graphs(args.input), args.output)
