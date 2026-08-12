#
#  makefile for wave equation project
#
# Compiler command
CC = gcc

#  compile/link flags, uncomment one, the top is for debugging
# CFLAGS = -g -Wall
CFLAGS = -O3 -Wall
#
# libraries: lapack, blas and standard math
#LIBS = -llapack -lblas -lgaussq -lm
#LFLAGS = -L/home/tausch/lib
LIBS = -llapack -lblas -lm

HEADERS = globals.h kernel.h mesh.h quadRule.h 
OBJ =  quadFFR.o polytope.o kernelWelle.o mesh.o quadRule.o solver.o

%.o: %.c $(HEADERS) Makefile
	$(CC) $(CFLAGS) -c -o $@ $<


welle: $(HEADERS) $(OBJ) welle.o Makefile 
	$(CC) $(CFLAGS) -o welle welle.o $(OBJ) $(LIBS) 


clean:
	rm *.o 


