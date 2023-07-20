default: compile run

compile:
	g++ -std=c++11 -I ./Libraries/include/eigen3 -I ./Libraries/include/eigen3/Eigen RSB/*.cpp Graph/*.cpp KL/*.cpp main.cpp -o main

run:
	./main

clean:
	rm -rf main