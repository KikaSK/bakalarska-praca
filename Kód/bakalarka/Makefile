CC = clang++
CCFLAGS = -std=c++17 -Wall -Wextra -g
GINAC = -lginac -lcln
BUILD_DIR = ./build
TEST_BUILD_DIR = ./test_build

all: run_main

reformat:
	clang-format -i -style=LlVM src/*.cpp src/*.h main.cpp

main: build/main.o build/function.o build/vector.o build/point.o build/triangle.o build/edge.o build/mesh.o
	$(CC) $(CCFLAGS) $(BUILD_DIR)/*.o -Isrc -o $(BUILD_DIR)/a.out $(GINAC)

run_main: main
	$(BUILD_DIR)/a.out

build/main.o: main.cpp src/basic_algorithm.h
	$(CC) $(CCFLAGS) -Isrc -c main.cpp -o $(BUILD_DIR)/main.o

build/function.o: src/function.cpp
	$(CC) $(CCFLAGS) -c src/function.cpp -o $(BUILD_DIR)/function.o

build/vector.o: src/vector.cpp
	$(CC) $(CCFLAGS) -c src/vector.cpp -o $(BUILD_DIR)/vector.o

build/point.o: src/point.cpp
	$(CC) $(CCFLAGS) -c src/point.cpp -o $(BUILD_DIR)/point.o

build/triangle.o: src/triangle.cpp
	$(CC) $(CCFLAGS) -c src/triangle.cpp -o $(BUILD_DIR)/triangle.o 

build/edge.o: src/edge.cpp
	$(CC) $(CCFLAGS) -c src/edge.cpp -o $(BUILD_DIR)/edge.o

build/mesh.o: src/mesh.cpp
	$(CC) $(CCFLAGS) -c src/mesh.cpp -o $(BUILD_DIR)/mesh.o


GOOGLETESDIR = ../../../../../Installs/googletest
FUSED_GTEST_DIR = $(GOOGLETESDIR)/googletest/src

# Paths to the fused gtest files.
FUSED_GTEST_H = $(FUSED_GTEST_DIR)/gtest/gtest.h
FUSED_GTEST_ALL_CC = $(FUSED_GTEST_DIR)/gtest/gtest-all.cc

# Where to find the sample test.
SAMPLE_DIR = ./test

# Where to find gtest_main.cc.
GTEST_MAIN_CC = $(GOOGLETESDIR)/googletest/src/gtest_main.cc

# Flags passed to the preprocessor.
# We have no idea here whether pthreads is available in the system, so
# disable its use.
CPPFLAGS += -I$(FUSED_GTEST_DIR) -DGTEST_HAS_PTHREAD=0

# Flags passed to the C++ compiler.
CXXFLAGS += -g -std=c++17

test_all : test1

test_check : test_all
	$(TEST_BUILD_DIR)/test1

test_clean :
	rm -rf  test1 *.o 
# $(FUSED_GTEST_DIR)

$(FUSED_GTEST_H) :
	$(FUSED_GTEST_DIR)/../scripts/fuse_gtest_files.py $(FUSED_GTEST_DIR)

$(FUSED_GTEST_ALL_CC):
	$(FUSED_GTEST_DIR)/../scripts/fuse_gtest_files.py $(FUSED_GTEST_DIR)

gtest-all.o : $(FUSED_GTEST_H) $(FUSED_GTEST_ALL_CC)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(FUSED_GTEST_DIR)/gtest/gtest-all.cc

gtest_main.o : $(FUSED_GTEST_H) $(GTEST_MAIN_CC)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(GTEST_MAIN_CC)

test1.o : $(SAMPLE_DIR)/test1.cc $(FUSED_GTEST_H)
	$(CXX) -no-pie $(CPPFLAGS) $(CXXFLAGS) -Isrc/ -c $(SAMPLE_DIR)/test1.cc $(GINAC)

test1 : test1.o gtest-all.o gtest_main.o build/point.o build/vector.o build/edge.o build/triangle.o build/function.o build/mesh.o
	$(CXX) -no-pie $(CPPFLAGS) $(CXXFLAGS) ./*.o ./build/point.o ./build/vector.o ./build/edge.o ./build/triangle.o ./build/function.o ./build/mesh.o -o $(TEST_BUILD_DIR)/$@ $(GINAC)

