CC = g++
CFLAGS = -g -Wall

OBJECTS = Basic_Ising.o

ISING: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm $(OBJECTS)
