CC = g++
CFLAGS = -g -Wall

OBJECTS = SFC_Ising.o sfc.o

ISING: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm $(OBJECTS)
