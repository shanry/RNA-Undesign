################################
# Makefile
#
# author: Tianshuo Zhou, Wei Yu Tang
# edited by: 09/2023
################################

CC=g++
DEPS=src/*.hpp src/LinearFold.h src/Utils/energy_parameter.h src/Utils/feature_weight.h src/Utils/intl11.h src/Utils/intl21.h src/Utils/intl22.h src/Utils/utility_v.h  src/Utils/utility_v_max.h src/Utils/utility.h 
CFLAGS=-std=c++11 -O3 -fopenmp
objects=bin/main bin/main4

all: main main_nosh

main: main.o utils.o comps.o eval.o
	$(CC)  $(CFLAGS)  bin/main.o bin/utils.o bin/comps.o bin/eval.o  -o bin/main

main_nosh: main_nosh.o utils.o comps.o eval_nosh.o # no special hp 
	$(CC)  $(CFLAGS)  bin/main_nosh.o bin/utils.o bin/comps.o bin/eval_nosh.o  -o bin/main_nosh

main_mac: main.o utils.o comps.o eval.o
	$(CC)  $(CFLAGS) -Wl,-ld_classic  bin/main.o bin/utils.o bin/comps.o bin/eval.o  -o bin/main

main_nosh_mac: main_nosh.o utils.o comps.o eval_nosh.o # no special hp
	$(CC)  $(CFLAGS) -Wl,-ld_classic  bin/main_nosh.o bin/utils.o bin/comps.o bin/eval_nosh.o  -o bin/main_nosh

main.o: src/main.cpp $(DEPS)
	$(CC) src/main.cpp -c $(CFLAGS) -Dlv -DSPECIAL_HP -Dis_cube_pruning -Dis_candidate_list -o bin/main.o

main_nosh.o: src/main.cpp $(DEPS)
	$(CC) src/main.cpp -c $(CFLAGS) -Dlv -Dis_cube_pruning -Dis_candidate_list -o bin/main_nosh.o

utils.o: src/utils.cpp src/utils.h
	$(CC) -c src/utils.cpp $(CFLAGS) -o bin/utils.o

comps.o: src/comps.cpp src/comps.h
	$(CC) -c src/comps.cpp $(CFLAGS) -o bin/comps.o

eval.o: src/eval.cpp src/eval.h
	$(CC) -c src/eval.cpp $(CFLAGS) -Dlv -DSPECIAL_HP -Dis_cube_pruning -Dis_candidate_list  -o bin/eval.o

eval_nosh.o: src/eval.cpp src/eval.h
	$(CC) -c src/eval.cpp $(CFLAGS) -Dlv -Dis_cube_pruning -Dis_candidate_list  -o bin/eval_nosh.o

clean:
	-rm $(objects)
