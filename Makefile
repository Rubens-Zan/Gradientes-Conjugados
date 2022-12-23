# Autor: Rubens Laszlo
# Data: 12/2022
# GRR 20206147 


CC = gcc
EXEC = cgSolver
CFLAG = -Wall -g -std=c99 -lm
MODULOS = resolverGradConj \
	sislin \
	utils 
 
OBJETOS = main.o $(addsuffix .o,$(MODULOS) $(CFLAG))

.PHONY: all clean

all: CGSOLVER

CGSOLVER: $(OBJETOS)
	$(CC) -o $(EXEC) $(OBJETOS) $(CFLAG)

clean:
	$(RM) $(OBJETOS) $(EXEC)