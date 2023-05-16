gcc_flags= -fopenmp -g -lm -pthread -O3 -Wall -Wextra -pedantic

EXE=gillespie

SRC=$(EXE).cpp

run: $(EXE)
	./$(EXE)
	python analyse_gil.py


$(EXE): $(SRC)
		g++ $(gcc_flags) $< -o $@

clean:
	rm -rf $(EXE)
