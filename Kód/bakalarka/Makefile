CC = clang++
CCFLAGS = -std=c++17 -Wall -Wextra -fPIC -O3 -march=native
GINAC = -lginac -lcln
BUILD_DIR = ./build
TEST_BUILD_DIR = ./test_build
SRC_DIR = ./src
SAMPLE_TEST_DIR = ./test

all: prepare run_main

clean:
	rm -rf $(BUILD_DIR) $(TEST_BUILD_DIR)

prepare:
	mkdir -p $(BUILD_DIR) $(TEST_BUILD_DIR)

reformat:
	clang-format -i -style=LlVM $(SRC_DIR)/*.cpp $(SRC_DIR)/*.h main.cpp $(SAMPLE_TEST_DIR)/*.cc

main: build/main.o build/algorithms.o build/basic_algorithm.o build/function.o build/bounding_box.o build/vector.o build/point.o build/triangle.o build/edge.o build/mesh.o
	$(CC) $(CCFLAGS) $(BUILD_DIR)/*.o -Isrc -o $(BUILD_DIR)/main.out $(GINAC)

run_main: main
	$(BUILD_DIR)/main.out

build/main.o: main.cpp
	$(CC) $(CCFLAGS) -Isrc -c main.cpp -o $(BUILD_DIR)/main.o

build/algorithms.o: $(SRC_DIR)/algorithms.cpp
	$(CC) $(CCFLAGS) -c $(SRC_DIR)/algorithms.cpp -o $(BUILD_DIR)/algorithms.o

build/basic_algorithm.o: $(SRC_DIR)/basic_algorithm.cpp
	$(CC) $(CCFLAGS) -c $(SRC_DIR)/basic_algorithm.cpp -o $(BUILD_DIR)/basic_algorithm.o

build/function.o: $(SRC_DIR)/function.cpp
	$(CC) $(CCFLAGS) -c $(SRC_DIR)/function.cpp -o $(BUILD_DIR)/function.o

build/vector.o: $(SRC_DIR)/vector.cpp
	$(CC) $(CCFLAGS) -c $(SRC_DIR)/vector.cpp -o $(BUILD_DIR)/vector.o

build/bounding_box.o: $(SRC_DIR)/bounding_box.cpp
	$(CC) $(CCFLAGS) -c $(SRC_DIR)/bounding_box.cpp -o $(BUILD_DIR)/bounding_box.o

build/point.o: $(SRC_DIR)/point.cpp
	$(CC) $(CCFLAGS) -c $(SRC_DIR)/point.cpp -o $(BUILD_DIR)/point.o

build/triangle.o: $(SRC_DIR)/triangle.cpp
	$(CC) $(CCFLAGS) -c $(SRC_DIR)/triangle.cpp -o $(BUILD_DIR)/triangle.o 

build/edge.o: $(SRC_DIR)/edge.cpp
	$(CC) $(CCFLAGS) -c $(SRC_DIR)/edge.cpp -o $(BUILD_DIR)/edge.o

build/mesh.o: $(SRC_DIR)/mesh.cpp
	$(CC) $(CCFLAGS) -c $(SRC_DIR)/mesh.cpp -o $(BUILD_DIR)/mesh.o


GOOGLETESDIR = ../../../../../Installs/googletest
FUSED_GTEST_DIR = $(GOOGLETESDIR)/googletest/src

# Paths to the fused gtest files.
FUSED_GTEST_H = $(FUSED_GTEST_DIR)/gtest/gtest.h
FUSED_GTEST_ALL_CC = $(FUSED_GTEST_DIR)/gtest/gtest-all.cc

# Where to find gtest_main.cc.
GTEST_MAIN_CC = $(GOOGLETESDIR)/googletest/src/gtest_main.cc

# Flags passed to the preprocessor.
# We have no idea here whether pthreads is available in the system, so
# disable its use.
CPPFLAGS += -I$(FUSED_GTEST_DIR) -DGTEST_HAS_PTHREAD=0

# Flags passed to the C++ compiler.
CXXFLAGS += -g -std=c++17

test : prepare test1

test_check : test
	$(TEST_BUILD_DIR)/test1

test_clean :
	rm -rf  test1 *.o 
# $(FUSED_GTEST_DIR)

$(FUSED_GTEST_H) :
	$(FUSED_GTEST_DIR)/../scripts/fuse_gtest_files.py $(FUSED_GTEST_DIR)

$(FUSED_GTEST_ALL_CC):
	$(FUSED_GTEST_DIR)/../scripts/fuse_gtest_files.py $(FUSED_GTEST_DIR)

$(TEST_BUILD_DIR)/gtest-all.o : $(FUSED_GTEST_H) $(FUSED_GTEST_ALL_CC)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(FUSED_GTEST_DIR)/gtest/gtest-all.cc -o $(TEST_BUILD_DIR)/gtest-all.o

$(TEST_BUILD_DIR)/gtest_main.o : $(FUSED_GTEST_H) $(GTEST_MAIN_CC)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(GTEST_MAIN_CC) -o $(TEST_BUILD_DIR)/gtest_main.o

$(TEST_BUILD_DIR)/test1.o : $(SAMPLE_TEST_DIR)/test1.cc $(FUSED_GTEST_H) src/algorithms.h src/basic_algorithm.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Isrc/ -c $(SAMPLE_TEST_DIR)/test1.cc $(GINAC) -o $(TEST_BUILD_DIR)/test1.o

test1 : $(TEST_BUILD_DIR)/test1.o $(TEST_BUILD_DIR)/gtest-all.o $(TEST_BUILD_DIR)/gtest_main.o build/point.o build/vector.o build/bounding_box.o build/edge.o build/triangle.o build/function.o build/mesh.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(TEST_BUILD_DIR)/gtest_main.o $(TEST_BUILD_DIR)/gtest-all.o $(TEST_BUILD_DIR)/test1.o ./build/point.o ./build/vector.o ./build/bounding_box.o ./build/edge.o ./build/triangle.o ./build/function.o ./build/mesh.o -o $(TEST_BUILD_DIR)/$@ $(GINAC)

