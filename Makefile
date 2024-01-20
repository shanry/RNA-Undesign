################################
# Makefile
#
# author: Tianshuo Zhou, Wei Yu Tang
# edited by: 09/2023
################################

CC=g++
DEPS=src/*.hpp  src/csv.cpp src/LinearFold.h src/Utils/energy_parameter.h src/Utils/feature_weight.h src/Utils/intl11.h src/Utils/intl21.h src/Utils/intl22.h src/Utils/utility_v.h  src/Utils/utility_v_max.h src/Utils/utility.h 
CFLAGS=-std=c++11 -O3 -fopenmp
objects=bin/main bin/main4

all: main main0

main: main.o utils.o comps.o eval.o
	$(CC)  $(CFLAGS) bin/main.o bin/utils.o bin/comps.o bin/eval.o  -o bin/main

main0: main0.o utils.o comps.o eval0.o # only special hp of triloop
	$(CC)  $(CFLAGS) bin/main0.o bin/utils.o bin/comps.o bin/eval0.o  -o bin/main0

main.o: src/main.cpp $(DEPS)
	$(CC) src/main.cpp -c $(CFLAGS) -Dlv -DSPECIAL_HP_4 -Dis_cube_pruning -Dis_candidate_list -o bin/main.o

main0.o: src/main.cpp $(DEPS)
	$(CC) src/main.cpp -c $(CFLAGS) -Dlv -Dis_cube_pruning -Dis_candidate_list -o bin/main0.o

utils.o: src/utils.cpp src/utils.h
	$(CC) -c src/utils.cpp $(CFLAGS) -o bin/utils.o

comps.o: src/comps.cpp src/comps.h
	$(CC) -c src/comps.cpp $(CFLAGS) -o bin/comps.o

eval.o: src/eval.cpp src/eval.h
	$(CC) -c src/eval.cpp $(CFLAGS) -Dlv -DSPECIAL_HP_4 -Dis_cube_pruning -Dis_candidate_list  -o bin/eval.o

eval0.o: src/eval.cpp src/eval.h
	$(CC) -c src/eval.cpp $(CFLAGS) -Dlv -Dis_cube_pruning -Dis_candidate_list  -o bin/eval0.o

clean:
	-rm $(objects)