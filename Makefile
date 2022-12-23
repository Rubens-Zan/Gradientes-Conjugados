# Autor: Rubens Laszlo
# Data: 12/2022
# GRR 20206147 


CC = gcc
EXEC = cgSolver
CFLAG = -Wall -g
MODULOS = resolvedorGradConjug \
	sislin \
	utils 
 
OBJETOS = main.o $(addsuffix .o,$(MODULOS))

.PHONY: all clean
lib_sislin.o: lib_sislin.c lib_sislin.h
	$(CC) $(CFLAG) -c lib_sislin.c -lm

resolvedorGradConjug: resolvedorGradConjug.c resolvedorGradConjug.h
	$(CC) $(CFLAG) -c resolvedorGradConjug.c -lm

sislin: sislin.c sislin.h
	$(CC) $(CFLAG) -c sislin.c -lm

utils: utils.c utils.h
	$(CC) $(CFLAG) -c utils.c -lm

all: CGSOLVER

CGSOLVER: $(OBJETOS)
	$(CC) -o $(EXEC) $(OBJETOS) $(CFLAG)

clean:
	$(RM) $(OBJETOS) $(EXEC)