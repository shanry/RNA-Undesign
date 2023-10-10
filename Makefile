################################
# Makefile
#
# author: Tianshuo Zhou, Wei Yu Tang
# edited by: 09/2023
################################

CC=g++
DEPS=src/Utils/energy_parameter.h src/Utils/feature_weight.h src/Utils/intl11.h src/Utils/intl21.h src/Utils/intl22.h src/Utils/utility_v.h src/Utils/utility.h 
CFLAGS=-std=c++11 -O3
objects=bin/eval

all: eval main

eval: src/eval.cpp $(DEPS) 
	mkdir -p bin
	$(CC) src/eval.cpp $(CFLAGS) -Dlv -o $(objects)

main: src/main.cpp #src/eval.cpp $(DEPS)
	$(CC) src/main.cpp src/eval.cpp $(CFLAGS) -o bin/main
clean:
	-rm $(objects)
