################################
# Makefile
#
# author: Tianshuo Zhou, Wei Yu Tang
# edited by: 09/2023
################################

CC=g++
DEPS=src/LinearFold.h src/Utils/energy_parameter.h src/Utils/feature_weight.h src/Utils/intl11.h src/Utils/intl21.h src/Utils/intl22.h src/Utils/utility_v.h  src/Utils/utility_v_max.h src/Utils/utility.h 
CFLAGS=-std=c++11 -O3 -fopenmp
objects=bin/main

all: main

main: src/main.cpp #src/eval.cpp $(DEPS)
	$(CC) src/main.cpp  $(CFLAGS) -Dlv -Dis_cube_pruning -Dis_candidate_list -o bin/main
clean:
	-rm $(objects)
