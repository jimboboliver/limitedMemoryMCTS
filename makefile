default: main.c
	gcc -o mcts main.c -lm
	gcc -o debug main.c -lm -g


debug: main.c

clean:
	del mcts.exe
	del debug.exe