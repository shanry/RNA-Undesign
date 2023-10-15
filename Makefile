################################
# Makefile
#
# author: Tianshuo Zhou, Wei Yu Tang
# edited by: 09/2023
################################

CC=g++
DEPS=src/fold.h src/LinearFoldEval.h src/LinearFold.h src/Utils/energy_parameter.h src/Utils/feature_weight.h src/Utils/intl11.h src/Utils/intl21.h src/Utils/intl22.h src/Utils/utility_v.h src/Utils/utility.h 
CFLAGS=-std=c++11 -O3
objects=bin/eval bin/main

all: main 

main: src/main.cpp   $(DEPS)
	$(CC) $(CFLAGS) -Dlv -Dis_cube_pruning -Dis_candidate_list src/main.cpp -o bin/main

clean:
	-rm $(objects)
