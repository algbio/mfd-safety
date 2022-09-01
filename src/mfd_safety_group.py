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


def solve_group_testing(graphs, group_tests, output_file):

    mfds = solve_instances(graphs)

    output = open(output_file, 'w+')

    for g, (mfd, group_test) in enumerate(zip(mfds, group_tests)):

        output.write(f"# graph {g}\n")

        try:
            model = build_ilp_model_avoiding_multiple_paths(mfd, len(mfd['solution']), group_test['paths'])
            model.optimize()

            for unsafe_path in get_unsafe_paths(model, group_test['paths']):
                output.write(' '.join(map(str, sorted(set([item for t in unsafe_path for item in t])))))
                output.write('\n')

        except gp.GurobiError as e:
            print(f'Error code {e.errno}: {e}', file=sys.stderr)

        except AttributeError:
            print('Encountered an attribute error', file=sys.stderr)

    output.close()


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
    requiredNamed.add_argument('-g', '--group-test', type=str, help='Input (group tests) filename', required=True)
    requiredNamed.add_argument('-o', '--output', type=str, help='Output filename', required=True)

    args = parser.parse_args()

    threads = args.threads
    if threads == 0:
        threads = os.cpu_count()
    print(f"INFO: Using {threads} threads for the Gurobi solver")

    solve_group_testing(*read_input(args.input, args.group_test), args.output)
