CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2 -g #-fsanitize=address
TARGET = projectA
TEST_TARGET = test
CWD = $(shell pwd)

# Directories
SRC_DIR = $(CWD)/src
TEST_DIR = $(CWD)/tests
BIN_DIR = $(CWD)/bin
OBJ_DIR = $(CWD)/obj
INC_DIR = $(CWD)/include
LIB_DIR = $(CWD)/lib
DEPS_DIR = $(CWD)/deps
ALGO_DIR = $(CWD)/algorithms

# Executables
EXE = $(BIN_DIR)/$(TARGET)
TEST_EXE = $(BIN_DIR)/$(TEST_TARGET)

# Library flags
LIB_FLAGS = -lz -lgssw -lm -lstdc++ -lgt_gwfa -lgwfa -ledlib #-fsanitize=address
LDFLAGS = -L$(LIB_DIR)

# LibrariesF
LIBS = $(LIB_DIR)/libgt_gwfa.a $(LIB_DIR)/libgssw.a

# Source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
TEST_SRCS = $(wildcard $(TEST_DIR)/*.cpp)
ALGO_SRCS = $(wildcard $(ALGO_DIR)/*.cpp)

# Object files
OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(filter-out $(SRC_DIR)/main.cpp, $(SRCS)))
MAIN_OBJ = $(OBJ_DIR)/main.o
TEST_OBJS = $(patsubst $(TEST_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(TEST_SRCS))
OBJS += $(patsubst $(ALGO_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(ALGO_SRCS))

# Header files
HEADERS = $(wildcard $(SRC_DIR)/*.hpp)
TEST_HEADERS = $(wildcard $(TEST_DIR)/*.hpp)
ALGO_HEADERS = $(wildcard $(ALGO_DIR)/*.hpp)

# Default target
all: $(EXE)

tests: $(TEST_EXE) run_tests

algorithms: $(LIBS)


# Linking step for the executable
$(EXE): $(LIBS) $(OBJS) $(MAIN_OBJ)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(LDFLAGS) -o $@ $(filter %.o, $^) -I$(INC_DIR) $(LIB_FLAGS)

$(TEST_EXE): $(LIBS) $(OBJS) $(TEST_OBJS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(LDFLAGS) -o $@ $(filter %.o, $^) -I$(INC_DIR) $(LIB_FLAGS)


# Compile src dir
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(INC_DIR)
	cp $(HEADERS) $(INC_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $< -I$(INC_DIR)

# Compile main file
$(OBJ_DIR)/main.o: $(SRC_DIR)/main.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $< -I$(INC_DIR)

# Compile test dir
$(OBJ_DIR)/%.o: $(TEST_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(INC_DIR)
	cp $(TEST_HEADERS) $(INC_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $< -I$(INC_DIR)

# Compile algorithm dir
$(OBJ_DIR)/%.o: $(ALGO_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(INC_DIR)
	@mkdir -p $(INC_DIR)/algorithms
	cp -r $(ALGO_HEADERS) $(INC_DIR)/algorithms
	$(CXX) $(CXXFLAGS) -c -o $@ $< -I$(INC_DIR)


# Only ever use either gt_gwfa or gssw.
# Create gt_gwfa library
$(LIB_DIR)/libgt_gwfa.a:
	@mkdir -p $(LIB_DIR)
	@mkdir -p $(INC_DIR)
	@mkdir -p $(INC_DIR)/gt_gwfa
	+ cd $(ALGO_DIR)/gt_gwfa && $(MAKE) && cp -r lib/* $(CWD)/lib && cp -r include/* $(INC_DIR)/gt_gwfa && cp -r src/*.h $(INC_DIR)/gt_gwfa
	+ rm -rf $(LIB_DIR)/libgssw.a

# Create gssw library
$(LIB_DIR)/libgssw.a:
	@mkdir -p $(LIB_DIR)
	@mkdir -p $(INC_DIR)
	@mkdir -p $(INC_DIR)/gssw
	+ cd $(ALGO_DIR)/gssw && $(MAKE) && cp -r lib/* $(CWD)/lib && cp -r src/*.h $(INC_DIR)/gssw && cp -r src/simde $(INC_DIR)/gssw


run_tests:
	$(shell $(TEST_EXE))

# Clean up
clean:
	@rm -rf $(OBJ_DIR) $(INC_DIR) $(BIN_DIR) $(LIB_DIR)
	cd $(ALGO_DIR)/gt_gwfa && $(MAKE) clean
	cd $(ALGO_DIR)/gssw && $(MAKE) clean

.PHONY: all clean copy_headers