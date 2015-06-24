# ======================
# Opcions d'optimització
# ======================
COMPILER=gcc
OPT=-g -Wall -std=c99
#OPT=-O3 -Wall -pedantic -std=c99 

# =========
# Utilitats
# =========
main_Lm : main_Lm.o mean.o
	$(COMPILER) -o main_Lm $(OPT) main_Lm.o mean.o -lm

# ===========
# Bibliotecas
# ===========
main_Lm.o : main_Lm.c
	$(COMPILER) -c $(OPT) main_Lm.c -lm
mean.o : mean.c
	$(COMPILER) -c $(OPT) mean.c -lm

# ======
# Neteja
# ======
clean :
	rm -f *.o
realclean : clean
	rm -f main_Lm


