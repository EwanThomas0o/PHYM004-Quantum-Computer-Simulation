CC = gcc
CFLAGS = -Wall -g -lgsl -lcblas 

all: test

quantcomp.o: quantcomp.c quantcomp.h
	$(CC) -c quantcomp.c
do:
	./quant
	
clean: 
	rm quant

test: test.c quantcomp.o
	$(CC) $(CFLAGS) -o $@ $^

