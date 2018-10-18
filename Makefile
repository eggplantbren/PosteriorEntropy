CXX = g++
FLAGS = -std=c++17 -O3 -march=native -DNDEBUG -Wall -Wextra -pedantic -I include

default:
	make -C InfoNest/cpp
	$(CXX) $(FLAGS) -c include/Demo.h
	$(CXX) $(FLAGS) -c include/Planet.h
	$(CXX) $(FLAGS) -c src/main.cpp
	$(CXX) -L InfoNest/cpp -o main main.o -linfonest
	rm -f *.o
