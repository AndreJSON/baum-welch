.DEFAULT_GOAL := a4
CC = g++
FLAGS = -std=c++11 -Wall -Wextra

a1:
	$(CC) $(FLAGS) -o a1.out a1.cpp
	@echo "----------DONE!----------"

a2:
	$(CC) $(FLAGS) -o a2.out a2.cpp
	@echo "----------DONE!----------"

a3:
	$(CC) $(FLAGS) -o a3.out a3.cpp
	@echo "----------DONE!----------"

a4:
	$(CC) $(FLAGS) -o a4.out a4.cpp
	@echo "----------DONE!----------"

clean:
	rm -f *.out *.o