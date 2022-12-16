# Autor: Rubens Laszlo
# Data: 12/2022
# GRR 20206147 


CC = gcc
EXEC = CG_SOLVER
CFLAG = -Wall -std=c99
MODULOS = opmatrizes \
	sislin \
	utils 
 
OBJETOS = main.o $(addsuffix .o,$(MODULOS))

.PHONY: all clean

all: CGSOLVER

CGSOLVER: $(OBJETOS)
	$(CC) -o $(EXEC) $(OBJETOS) $(ALLEGRO) $(CFLAG)

clean:
	$(RM) $(OBJETOS)

purge:
	$(RM) $(OBJETOS) $(EXEC)