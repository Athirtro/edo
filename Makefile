CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -O2
LDLIBS = -lm

TESTS = artifacts/test_model artifacts/test_parameters artifacts/test_equilibria artifacts/test_collapse artifacts/test_sensitivity

.PHONY: test reference clean

test: $(TESTS)
	./artifacts/test_model
	./artifacts/test_parameters
	./artifacts/test_equilibria
	./artifacts/test_collapse
	./artifacts/test_sensitivity

artifacts/test_model: artifacts/test_model.c src/model.c src/model.h
	$(CC) $(CFLAGS) -Isrc artifacts/test_model.c src/model.c -o $@ $(LDLIBS)

artifacts/test_parameters: artifacts/test_parameters.c src/parameters.c src/parameters.h src/model.c src/model.h
	$(CC) $(CFLAGS) -Isrc artifacts/test_parameters.c src/parameters.c src/model.c -o $@ $(LDLIBS)

artifacts/test_equilibria: artifacts/test_equilibria.c src/analysis.c src/analysis.h src/parameters.c src/parameters.h src/model.c src/model.h
	$(CC) $(CFLAGS) -Isrc artifacts/test_equilibria.c src/analysis.c src/parameters.c src/model.c -o $@ $(LDLIBS)

artifacts/test_collapse: artifacts/test_collapse.c src/collapse.c src/collapse.h src/analysis.c src/analysis.h src/parameters.c src/parameters.h src/model.c src/model.h
	$(CC) $(CFLAGS) -Isrc artifacts/test_collapse.c src/collapse.c src/analysis.c src/parameters.c src/model.c -o $@ $(LDLIBS)

artifacts/test_sensitivity: artifacts/test_sensitivity.c src/sensitivity.c src/sensitivity.h src/collapse.c src/collapse.h src/analysis.c src/analysis.h src/parameters.c src/parameters.h src/model.c src/model.h
	$(CC) $(CFLAGS) -Isrc artifacts/test_sensitivity.c src/sensitivity.c src/collapse.c src/analysis.c src/parameters.c src/model.c -o $@ $(LDLIBS)

artifacts/reference_rk4: artifacts/reference_rk4.c src/collapse.c src/collapse.h src/parameters.c src/parameters.h src/model.c src/model.h
	$(CC) $(CFLAGS) -Isrc artifacts/reference_rk4.c src/collapse.c src/analysis.c src/parameters.c src/model.c -o $@ $(LDLIBS)

reference: artifacts/reference_rk4
	./artifacts/reference_rk4

clean:
	rm -f $(TESTS) artifacts/reference_rk4
