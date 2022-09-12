#!/usr/bin/env python
# coding: utf-8

import os
import argparse
import time
from mfd_safety import read_input, solve_instances, is_safe


def write_graph(g, mfd, output):
    n_nodes = len(mfd["vertices"])
    output.write(f"# graph {g}\n{n_nodes}\n")
    for edge in mfd["graph"].edges():
        weight = mfd["graph"].edges[edge[0], edge[1]]["weight"]
        output.write(f"{edge[0]} {edge[1]} {weight}\n")


def write_instances(graphs, output_file):
    mfds = solve_instances(graphs)

    output = open(output_file, 'w+')

    counter = 0

    for g, mfd in enumerate(mfds):
        counter += 1
        if counter % 100 == 0:
            print(counter)
        paths = mfd['solution']
        # first, let's check whether the mfd we got already is safe.
        mfd_unique = True
        for path in paths:
            if not is_safe(mfd, len(paths), path):
                mfd_unique = False
                break

        if not mfd_unique:
            # write
            write_graph(g, mfd, output)

    output.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="""
        Finds instances with a non-unique set of MFD paths.
        This script uses the Gurobi ILP solver.
        """,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-t', '--threads', type=int, default=0,
                        help='Number of threads to use for the Gurobi solver; use 0 for all threads (default 0).')

    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-i', '--input', type=str, help='Input filename', required=True)
    requiredNamed.add_argument('-o', '--output', type=str, help='Output filename', required=True)

    args = parser.parse_args()

    threads = args.threads
    if threads == 0:
        threads = os.cpu_count()
    print(f"INFO: Using {threads} threads for the Gurobi solver")

    start_time = time.time()
    write_instances(read_input(args.input), args.output)
    print(f"wall clock time: {time.time() - start_time}")
