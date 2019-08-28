CC = g++
CFLAGS = -g -Wall

OBJECTS = SFC_Ising.o sfc.o

ISING: $(OBJECTS)
	$(CC) -std=c++11 $(CFLAGS) -o $@ $^ -lm

clean:
	rm $(OBJECTS)
