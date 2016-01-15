# ----------------------------------
# FLAGS AND COMPILER TO USE IN LINUX
# ----------------------------------


CC=gcc-5
CCLLVM=gcc
CFLAGS=-c -g -DCLI -std=c99 -Wall -pedantic  
LFLAGS= -lm -L/usr/local/lib 

EXE=isr3
SRC=covarTree.c sweepTree.c  
OBJ=covarTree.o sweepTree.o 


$(EXE): $(OBJ) 
	$(CC) $(OBJ) $(LFLAGS) -o $(EXE) 


$(OBJ): $(SRC) 
#	$(CCLLVM) $(CFLAGS) $(SRC) 
#	rm -f $(OBJ) covarTree.o sweepTree.o 
	$(CC) $(CFLAGS) $(SRC) 

clean:
	rm -f $(EXE) $(OBJ) isr3 covarTree.o sweepTree.o 
