# Autor: Rubens Laszlo
# Data: 12/2022
# GRR 20206147 
CC=gcc
CFLAGS=-Wall -g
SISLIN=-DSISLIN
PROG=cgSolver
OBJS=resolvedorGradConjug.o sislin.o utils.o $(PROG).o

all: $(PROG)

sislin: CFLAGS += $(SISLIN)
sislin: all

resolvedorGradConjug.o: resolvedorGradConjug.c resolvedorGradConjug.h utils.h
	$(CC) $(CFLAGS) -c resolvedorGradConjug.c -lm

sislin.o: sislin.c sislin.h resolvedorGradConjug.h utils.h
	$(CC) $(CFLAGS) -c sislin.c -lm

utils.o: utils.c utils.h
	$(CC) $(CFLAGS) -c utils.c -lm

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm -f *.o *.bak $(PROG)
