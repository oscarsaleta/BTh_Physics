# ======================
# Opcions d'optimització
# ======================
COMPILER=gcc
OPT=-g -Wall -std=c99
#OPT=-O3 -Wall -pedantic -std=c99 

# =========
# Utilitats
# =========
main_L : main_L.o momang.o
	$(COMPILER) -o main_L $(OPT) main_L.o momang.o -lm

# ===========
# Bibliotecas
# ===========
main_L.o : main_L.c
	$(COMPILER) -c $(OPT) main_L.c -lm
momang.o : momang.c
	$(COMPILER) -c $(OPT) momang.c -lm

# ======
# Neteja
# ======
clean :
	rm -f *.o
realclean : clean
	rm -f main_L 


