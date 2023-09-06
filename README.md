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

2. Compile code

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

- MLRSB: Sequential MLRSB
- pMLRSB: Sequential pMLRSB
- Parallel_MLRSB: Parallel MLRSB
- Parallel_pMLRSB: Parallel pMLRSB
- RSB: Sequential RSB
- pRSB: Sequential pRSB
- Parallel_RSB: Parallel RSB
- Parallel_pRSB: Parallel pRSB

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

- MLRSB.c
- Parallel_MLRSB.c
- RSB.c
- Parallel_RSB.c
- MLRSB.h
- Parallel_MLRSB.h
- RSB.h
- Parallel_RSB.h

### utils

It contains all utility functions.

- utils.c
- utils.h

### Graph

It contains both functions related to graphs.

- Graph.c
- Graph.h

<!-- Files' Content -->

## Files' Content

#### MLRSB.c

TODO

#### Parallel_MLRSB.c

TODO

#### RSB.c

TODO

#### Parallel_RSB.c

TODO

#### utils.c

TODO

#### Graph.c

TODO
