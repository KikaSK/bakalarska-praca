# https://devhints.io/makefile
CC = g++
CCFLAGS = -std=c++17 -Wall -Wextra -fPIC -O2 #-O3 -march=native
GINAC = -lginac -lcln
BUILD_DIR = build
TEST_BUILD_DIR = test_build
SRC_DIR = src
SAMPLE_TEST_DIR = test
OUT_DIR = outputs
dependencies := $(wildcard src/*.cpp)
build_dependencies := $(dependencies:src/%.cpp=build/%.o)

all: run_main

clean:
	rm -rf $(BUILD_DIR) $(TEST_BUILD_DIR)

prepare:
	mkdir -p $(BUILD_DIR) $(TEST_BUILD_DIR) $(OUT_DIR)

reformat:
	clang-format -i -style=LlVM $(SRC_DIR)/*.cpp $(SRC_DIR)/*.h main.cpp measure.cpp $(SAMPLE_TEST_DIR)/*.cc

run_main: clean prepare main
	$(BUILD_DIR)/main.out

run_measure: clean prepare measure
	$(BUILD_DIR)/measure.out

run_test: clean prepare test
	$(BUILD_DIR)/test.out

main: build/main.o $(build_dependencies)
	$(CC) $(CCFLAGS) $^ -I$(SRC_DIR) -o $(BUILD_DIR)/main.out $(GINAC)

measure: build/measure.o $(build_dependencies)
	$(CC) $(CCFLAGS) $^ -I$(SRC_DIR) -o $(BUILD_DIR)/measure.out $(GINAC)

test: build/test.o $(build_dependencies)
	$(CC) $(CCFLAGS) $^ -I$(SRC_DIR) -o $(BUILD_DIR)/test.out $(GINAC)

build/%.o: %.cpp
	$(CC) $(CCFLAGS) -I$(SRC_DIR) -c $^ -o $@

build/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CCFLAGS) -c $< -o $@