gcc_flags= -fopenmp -g -lm -pthread -O3 -Wall -Wextra -pedantic

EXE=bruss

SRC=$(EXE).cpp

run: $(EXE)
	./$(EXE)


$(EXE): $(SRC)
		g++ $(gcc_flags) $< -o $@

clean:
	rm -rf $(EXE)
