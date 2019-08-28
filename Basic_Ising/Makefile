CC = g++
CFLAGS = -g -Wall

OBJECTS = Basic_Ising.o

ISING: $(OBJECTS)
	$(CC) $(CFLAGS) -std=c++11 -o $@ $^ -lm

clean:
	rm $(OBJECTS)
