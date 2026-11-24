CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -O2
LDLIBS = -lm

.PHONY: test clean

test: artifacts/test_model
	./artifacts/test_model

artifacts/test_model: artifacts/test_model.c src/model.c src/model.h
	$(CC) $(CFLAGS) -Isrc artifacts/test_model.c src/model.c -o $@ $(LDLIBS)

clean:
	rm -f artifacts/test_model
