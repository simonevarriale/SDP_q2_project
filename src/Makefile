default: compile run

compile_graph_generator:
	g++ -std=c++20 Generator/graph_generator.cpp -o graph_generator

run:
	./graph_partitioner

run_args:
	./graph_partitioner $(file) $(algorithm) $(partitions) $(threads)

compile:
	g++ -std=c++20 -I ./Libraries/include/eigen3 RSB/*.cpp Graph/*.cpp main.cpp utils/*.cpp -o graph_partitioner

clean:
	rm -rf graph_partitioner

generate_graph:
	./graph_generator $(n) $(e)

get_data:
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' RSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' MLRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' pRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' pRSB 4 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' pRSB 8 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' pMLRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' pMLRSB 4 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' pMLRSB 8 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_RSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_RSB 2 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_RSB 2 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_MLRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_MLRSB 2 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_MLRSB 2 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_pRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_pRSB 2 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_pRSB 2 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_pRSB 4 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_pRSB 4 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_pRSB 4 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_pRSB 8 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_pRSB 8 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_pRSB 8 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_pMLRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_pMLRSB 2 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_pMLRSB 2 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_pMLRSB 4 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_pMLRSB 4 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_pMLRSB 4 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_pMLRSB 8 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_pMLRSB 8 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_50_128.txt' Parallel_pMLRSB 8 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' RSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' MLRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' pRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' pRSB 4 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' pRSB 8 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' pMLRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' pMLRSB 4 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' pMLRSB 8 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_RSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_RSB 2 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_RSB 2 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_MLRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_MLRSB 2 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_MLRSB 2 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_pRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_pRSB 2 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_pRSB 2 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_pRSB 4 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_pRSB 4 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_pRSB 4 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_pRSB 8 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_pRSB 8 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_pRSB 8 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_pMLRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_pMLRSB 2 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_pMLRSB 2 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_pMLRSB 4 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_pMLRSB 4 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_pMLRSB 4 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_pMLRSB 8 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_pMLRSB 8 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_100_256.txt' Parallel_pMLRSB 8 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' RSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' MLRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' pRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' pRSB 4 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' pRSB 8 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' pMLRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' pMLRSB 4 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' pMLRSB 8 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_RSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_RSB 2 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_RSB 2 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_MLRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_MLRSB 2 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_MLRSB 2 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_pRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_pRSB 2 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_pRSB 2 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_pRSB 4 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_pRSB 4 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_pRSB 4 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_pRSB 8 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_pRSB 8 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_pRSB 8 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_pMLRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_pMLRSB 2 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_pMLRSB 2 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_pMLRSB 4 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_pMLRSB 4 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_pMLRSB 4 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_pMLRSB 8 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_pMLRSB 8 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_250_640.txt' Parallel_pMLRSB 8 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' RSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' MLRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' pRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' pRSB 4 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' pRSB 8 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' pMLRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' pMLRSB 4 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' pMLRSB 8 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_RSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_RSB 2 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_RSB 2 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_MLRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_MLRSB 2 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_MLRSB 2 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_pRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_pRSB 2 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_pRSB 2 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_pRSB 4 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_pRSB 4 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_pRSB 4 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_pRSB 8 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_pRSB 8 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_pRSB 8 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_pMLRSB 2 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_pMLRSB 2 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_pMLRSB 2 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_pMLRSB 4 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_pMLRSB 4 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_pMLRSB 4 8
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_pMLRSB 8 2
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_pMLRSB 8 4
	./graph_partitioner '/Users/gerardomaruotti/Developer/SDP_q2_project/src/data/graph_500_1024.txt' Parallel_pMLRSB 8 8