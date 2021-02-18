all:
	g++ -std=c++17 src/*.cpp main.cpp -Isrc -o a.out -lginac -lcln && time ./a.out

reformat:
	clang-format -i -style=LlVM src/*.cpp src/*.h
