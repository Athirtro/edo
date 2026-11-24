CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -O2
LDLIBS = -lm

TESTS = artifacts/test_model artifacts/test_parameters

.PHONY: test clean

test: $(TESTS)
	./artifacts/test_model
	./artifacts/test_parameters

artifacts/test_model: artifacts/test_model.c src/model.c src/model.h
	$(CC) $(CFLAGS) -Isrc artifacts/test_model.c src/model.c -o $@ $(LDLIBS)

artifacts/test_parameters: artifacts/test_parameters.c src/parameters.c src/parameters.h src/model.c src/model.h
	$(CC) $(CFLAGS) -Isrc artifacts/test_parameters.c src/parameters.c src/model.c -o $@ $(LDLIBS)

clean:
	rm -f $(TESTS)
