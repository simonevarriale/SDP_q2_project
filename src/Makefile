default: compile_eigen run clean

compile:
	g++ -std=c++20 Graph/*.cpp KL/*.cpp main.cpp -o graph_partitioner

run:
	./graph_partitioner

run_args:
	./graph_partitioner $(file) $(algorithm) $(partitions) $(threads)

compile_eigen:
	g++ -std=c++20 -I ./Libraries/include/eigen3 RSB/*.cpp Graph/*.cpp KL/*.cpp main.cpp utils/*.cpp -o graph_partitioner

clean:
	rm -rf graph_partitioner

generate_graph:
	g++ -std=c++20 Generator/graph_generator.cpp -o graph_generator
	./graph_generator $(n) $(e)
	rm -rf ./Generator/graph_generator