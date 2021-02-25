CC = g++
CCFLAGS = -std=c++17 -Wall -Wextra -g
GINAC = -lginac -lcln

all: run_main

reformat:
	clang-format -i -style=LlVM src/*.cpp src/*.h main.cpp

main: main.o function.o vector.o point.o triangle.o edge.o mesh.o
	$(CC) $(CCFLAGS) ./build/*.o -Isrc -o ./build/a.out $(GINAC)

run_main: main
	./build/a.out

main.o: main.cpp
	$(CC) $(CCFLAGS) -Isrc -c main.cpp -o build/main.o

function.o: src/function.cpp
	$(CC) $(CCFLAGS) -c src/function.cpp -o build/function.o

vector.o: src/vector.cpp
	$(CC) $(CCFLAGS) -c src/vector.cpp -o build/vector.o

point.o: src/point.cpp
	$(CC) $(CCFLAGS) -c src/point.cpp -o build/point.o

triangle.o: src/triangle.cpp
	$(CC) $(CCFLAGS) -c src/triangle.cpp -o build/triangle.o 

edge.o: src/edge.cpp
	$(CC) $(CCFLAGS) -c src/edge.cpp -o build/edge.o

mesh.o: src/mesh.cpp
	$(CC) $(CCFLAGS) -c src/mesh.cpp -o build/mesh.o
