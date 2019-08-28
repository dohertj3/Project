CC = mpiCC

CFLAGS = -Wall

OBJECTS = check_ising.o

ising: $(OBJECTS)
	$(CC) -std=c++11 $(CFLAGS) -o $@ $^ -lm

clean:
	rm $(OBJECTS)
