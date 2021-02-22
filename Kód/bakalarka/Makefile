CC = g++
CCFLAGS = -std=c++17
GINAC = -lginac -lcln

all: main

main: main.o function.o vector.o point.o triangle.o edge.o mesh.o
	$(CC) $(CCFLAGS) build/*.o -Isrc -o a.out $(GINAC) && ./a.out

main.o: main.cpp
	$(CC) $(CCFLAGS) -Isrc -c main.cpp -o build/main.o $(GINAC)

function.o: src/function.cpp
	$(CC) $(CCFLAGS) -c src/function.cpp -o build/function.o $(GINAC)

vector.o: src/vector.cpp
	$(CC) $(CCFLAGS) -c src/vector.cpp -o build/vector.o $(GINAC)

point.o: src/point.cpp
	$(CC) $(CCFLAGS) -c src/point.cpp -o build/point.o $(GINAC)

triangle.o: src/triangle.cpp
	$(CC) $(CCFLAGS) -c src/triangle.cpp -o build/triangle.o $(GINAC)

edge.o: src/edge.cpp
	$(CC) $(CCFLAGS) -c src/edge.cpp -o build/edge.o $(GINAC)

mesh.o: src/mesh.cpp
	$(CC) $(CCFLAGS) -c src/mesh.cpp -o build/mesh.o $(GINAC)


reformat:
	clang-format -i -style=LlVM src/*.cpp src/*.h main.cpp
