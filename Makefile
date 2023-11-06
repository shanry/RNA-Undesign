################################
# Makefile
#
# author: Tianshuo Zhou, Wei Yu Tang
# edited by: 09/2023
################################

CC=g++
DEPS=src/LinearFold.h src/Utils/energy_parameter.h src/Utils/feature_weight.h src/Utils/intl11.h src/Utils/intl21.h src/Utils/intl22.h src/Utils/utility_v.h  src/Utils/utility_v_max.h src/Utils/utility.h 
CFLAGS=-std=c++11 -O3 -fopenmp
objects=bin/main bin/main4

all: main main0

main: src/main.cpp $(DEPS)
	$(CC) src/main.cpp  $(CFLAGS) -Dlv -DSPECIAL_HP_4 -Dis_cube_pruning -Dis_candidate_list -o bin/main

main0: src/main.cpp $(DEPS) # only special hp of triloop
	$(CC) src/main.cpp  $(CFLAGS) -Dlv -Dis_cube_pruning -Dis_candidate_list -o bin/main0

clean:
	-rm $(objects)
