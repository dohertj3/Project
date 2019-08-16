CC = gcc

CFLAGS = -g -Wall

OBJECTS = test.o hilbert.o

TEST: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm $(OBJECTS)
