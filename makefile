CC = gcc

quant: quantcomp.c
	$(CC) -o quant quantcomp.c -lgsl -lcblas -Wall
do:
	./quant
	
clean: 
	rm quant