# ----------------------------------
# FLAGS AND COMPILER TO USE IN LINUX
# ----------------------------------


CC=gcc-5
CFLAGS=-c -g -DTEST_REBUILDCOVAR -std=c99 -Wall -pedantic  
#CFLAGS=-c -g -DTEST_SWEEPTREE -std=c99 -Wall -pedantic  
LFLAGS= -lm -L/usr/local/lib 

EXE=main
SRC= rebuildCovar.c sweepTree.c covarTree.c 
OBJ= rebuildCovar.o sweepTree.o covarTree.o 



$(EXE): $(OBJ) 
	$(CC) $(OBJ) $(LFLAGS) -o $(EXE) 

$(OBJ): $(SRC) 
	$(CC) $(CFLAGS) $(SRC) 

clean:
	rm -f $(EXE) $(OBJ) 

