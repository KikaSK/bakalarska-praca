all:
	g++ src/*.cpp main.cpp -Isrc -o a.out -lginac -lcln && ./a.out

reformat:
	clang-format -i -style=LlVM src/*.cpp src/*.h
