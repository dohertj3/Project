CC = mpicc

CFLAGS = -Wall

OBJECTS = check_ising.o

ising: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm $(OBJECTS)
