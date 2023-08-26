default: compile run clean

compile:
	g++ -std=c++20 Graph/*.cpp KL/*.cpp main.cpp -o main

run:
	./main

compile_eigen:
	g++ -std=c++20 -I ./Libraries/include/eigen3 -I ./Libraries/include/eigen3/Eigen RSB/*.cpp Graph/*.cpp KL/*.cpp main.cpp -o main

clean:
	rm -rf main