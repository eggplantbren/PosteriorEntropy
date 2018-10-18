CXX = g++
FLAGS = -std=c++17 -O2 -g -Wall -Wextra -pedantic

default:
	make -C InfoNest/cpp
	$(CXX) $(FLAGS) -c include/Demo.h
	$(CXX) $(FLAGS) -c include/Planet.h
