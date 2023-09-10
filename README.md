# Graph Partitioning

<!-- ABOUT THE PROJECT -->

## About The Project

In computer science, it is common to work with large graphs that need to be divided into non-overlapping sections, known as partitions. This process is called p-way partitioning. Our aim is to minimize the sum of the weights of edges that cross between partitions and balance the sum of weights of nodes in each partition.

<!-- GETTING STARTED -->

## Getting Started

The commands to run the entire project are in `Makefile` file.

### Installation

There are three fundamental commands:

1. Delete previous object files

```sh
make clean
```

2. Compile code (using g++ compiler)

```sh
make compile
```

3. Launch program

```sh
./graph_partitioner input_file algorithm_name partitions_number threads_number
```

4. Create graph file

```sh
./graph_generator nodes_number edges_number
```

This command will create a graph file in the data folder with the name graph_nodes_number_edges_number.txt. The weights of nodes and edges are random numbers between 1 and 10.

#### Input file

- ./data/graph_50_128.txt
- ./data/graph_100_256.txt
- ./data/graph_250_640.txt
- ./data/graph_500_1024.txt

#### Algorithm name

- MLRSB: Sequential Multi-Level Recursive Spectral Bisection
- pMLRSB: Sequential p-way Multi-Level Recursive Spectral Bisection
- Parallel_MLRSB: Parallel Multi-Level Recursive Spectral Bisection
- Parallel_pMLRSB: Parallel p-way Multi-Level Recursive Spectral Bisection
- RSB: Sequential Recursive Spectral Bisection
- pRSB: Sequential p-way Recursive Spectral Bisection
- Parallel_RSB: Parallel Recursive Spectral Bisection
- Parallel_pRSB: Parallel p-way Recursive Spectral Bisection

## Usage

```sh
 ./graph_partitioner  ./data/graph_50_128.txt MLRSB # Sequential MLRSB

./graph_partitioner  ./data/graph_50_128.txt pMLRSB 4 # pMLRSB with 4 partitions

./graph_partitioner  ./data/graph_50_128.txt Parallel_pMLRSB 4 4 # Parallel pMLRSB with 4 partitions and 4 threads

./graph_generator 50 128 # Create a graph file with 50 nodes and 128 edges
```

<!-- Folders' Content -->

## Folders' Content

There are three main folders.

### RSB

It contains all the functions related to RSB and MLRSB algorithm.

- MLRSB.cpp
- Parallel_MLRSB.cpp
- RSB.cpp
- Parallel_RSB.cpp
- MLRSB.h
- Parallel_MLRSB.h
- RSB.h
- Parallel_RSB.h

### utils

It contains all utility functions.

- utils.cpp
- utils.h

### Graph

It contains both functions related to graphs.

- Graph.cpp
- Graph.h

<!-- Files' Content -->

## Files' Content

#### MLRSB.cpp

It contains the main function `MLRSB` that implements the MLRSB algorithm. It takes as input the graph and it returns the best partitioning found. It also contains the function `pMLRSB` that implements the pMLRSB algorithm. It takes as input the graph and the number of partitions and it returns the best partitioning found.

#### Parallel_MLRSB.cpp

It contains the main function `Parallel_MLRSB` that implements the Parallel MLRSB algorithm. It takes as input the graph and the number of threads and it returns the best partitioning found. It also contains the function `Parallel_pMLRSB` that implements the Parallel pMLRSB algorithm. It takes as input the graph, the number of partitions and the number of threads and it returns the best partitioning found.

#### RSB.cpp

It contains the main function `RSB` that implements the RSB algorithm. It takes as input the graph and it returns the best partitioning found. It also contains the function `pRSB` that implements the pRSB algorithm. It takes as input the graph and the number of partitions and it returns the best partitioning found.

#### Parallel_RSB.cpp

It contains the main function `Parallel_RSB` that implements the Parallel RSB algorithm. It takes as input the graph and the number of threads and it returns the best partitioning found. It also contains the function `Parallel_pRSB` that implements the Parallel pRSB algorithm. It takes as input the graph, the number of partitions and the number of threads and it returns the best partitioning found.

#### utils.cpp

It contains a set of utility functions. The most important are: `read_graph` that reads the graph from a text file, `uncoarsening` that performs the uncoarsening phase of the algorithm and other functions to calculate the quality of the partitioning.

#### Graph.cpp

It contains the class `Graph` that represents the graph. It contins functions to calculate the adjacency matrix, degree matrix and other utility functions to add/remove nodes and edges. Finally, it contains some functions to help with the coarsening phase.
