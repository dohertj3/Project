CC = mpiCC

CFLAGS = -Wall

OBJECTS = ParSFC_Ising.o ParSFC.o

ising: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ -lm -std=c++11

clean:
	rm $(OBJECTS)
