# mfd-safety

mfd-safety is a tool reporting maximal safe paths for minimum flow decompositions (mfd) using ILP calls,
and implementing several optimization to reduce the number of ILP calls or their size (number of variables/constrains).

# Installation

- Install [miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
- Clone our this repository and `cd` to the corresponding folder
- `conda env create -f conda_environment.yml`

# Run

To run the project activate the conda environment you created during installation (`conda activate mfd-safety`) and use
`python` to execute the `mfd_safety.py` file.

As an example you can try:

`python ./src/mfd_safety.py -i ./example_inputs/example.graph -o ./example_inputs/example.safe`

## Input

- The input is a file containing a sequence of (directed) acyclic flow graphs separated by lines starting with `#`.
- The first line of each flow graph contains the number of vertices of the graph, after this every flow edge from vertex
`u` to  `v` carrying `f` flow is represented in a separated line in the format `u v f`.
- Vertices must be integers following a topological order of the graph.
- An example of such a format can be found in `./example_inputs/example.graph`.
- To specify subpath constraints, after the specification of each graph add a line 'subpaths' followed by the subpaths, one on each line. A subpath is specified as the corresponding sequence of vertices (separated by spaces) followed by '1.0', an example can be found in `./example_inputs/example_subpaths.graph`

## Output

- The output is a file containing a sequence of maximal safe paths separated by lines starting with `#` (one per flow
graph in the input).
- Each line contains a maximal safe path as the corresponding sequence of vertices (starting with -1).
- An example of such a format can be found in `./example_inputs/example.safe`.
- The maximality property of the safe paths is only guaranteed within a path in a flow decomposition, to obtain global
maximal safe path one must remove subpaths, for our experiments we used [this code](https://github.com/algbio/flow-decomposition-safety/blob/main/src/cpp-scripts/acTrie.cpp).

## Parameters

- `-i <path to input file>`. Mandatory.
- `-o <path to locate output>`. Mandatory.
- `-stats` Output stats to file <output>.stats
- `-t <n>` Use n threads for the Gurobi solver; use 0 for all threads (default 0).
- `-ilptb <n>` Maximum time (in seconds) that the ilp solver is allowed to take when computing safe paths for one flow graph.
If the solver takes more than n seconds, then safe for (all) flow decompositions is reported instead.
- `-uef` Uses excess flow to save ILP calls.
- `-uy2v` Use Y2V contraction on the flow graphs to reduce the ILP size.
- `-s/es/rs/ess/esl/rss/rsl {scan, bin_search, exp_search, rep_exp_search}` When running the two-finger algorithm applied
the specified strategy to extend/reduce the current safe interval.
- `-st/est/rst <n>` When running the two-finger algorithm run the `small strategy` when the search space is less than n
and the `large strategy` otherwise.
- `-ugtd/-ugbu` Run a group testing algorithm (top down or bottom up) instead of two-finger.


# Our experiments

## Commands

### Two-Finger

`python ./src/mfd_safety.py -i <input> -o <output> -t 12 -uy2v -uef [-ilptb 120]`

### Two-Finger Optimized Search

`python ./src/mfd_safety.py -i <input> -o <output> -t 12 -uy2v -uef -st 8 -esl bin_search -ess exp_search -rsl bin_search -rss exp_search [-ilptb 120]`

### Group Testing Bottom-Up

`python ./src/mfd_safety.py -i <input> -o <output> -t 12 -uy2v -uef -ugbu [-ilptb 120]`

### Group Testing Top-Down

`python ./src/mfd_safety.py -i <input> -o <output> -t 12 -uy2v -uef -ugtd [-ilptb 120]`

## Datasets

The datasets can be found in Zenodo at: [https://doi.org/10.5281/zenodo.8275700](https://doi.org/10.5281/zenodo.8275700
)
