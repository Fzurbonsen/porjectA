CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2 -g -fsanitize=address
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
LIB_FLAGS = -lz -lgssw -lm -lstdc++ -lgwfa -ledlib -lksw2 -lcsswl -lgnwa -fsanitize=address # -labpoa -lvargas_lib
LDFLAGS = -L$(LIB_DIR)

# Libraries
LIBS = $(LIB_DIR)/libgssw.a
LIBS += $(LIB_DIR)/libgwfa.a
LIBS += $(LIB_DIR)/libksw2.a
LIBS += $(LIB_DIR)/libcsswl.a
LIBS += $(LIB_DIR)/libedlib.a
# LIBS += $(LIB_DIR)/libabpoa.a
# LIBS += $(LIB_DIR)/libvargas_lib.a
LIBS += $(LIB_DIR)/libgnwa.a

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

# Ensure headers are prepared
prepare_headers:
	@mkdir -p $(INC_DIR)/nlohmann
	+ cd $(DEPS_DIR)/json/include && cp -r ./nlohmann/* $(INC_DIR)/nlohmann
	@mkdir -p $(INC_DIR)
	@mkdir -p $(INC_DIR)/algorithms
	@mkdir -p $(INC_DIR)/gssw
	@mkdir -p $(INC_DIR)/gwfa
	@mkdir -p $(INC_DIR)/ksw2
	@mkdir -p $(INC_DIR)/csswl
	@mkdir -p $(INC_DIR)/edlib
	@mkdir -p $(INC_DIR)/abPOA
	@mkdir -p $(INC_DIR)/vargas
	@mkdir -p $(INC_DIR)/vargas/htslib
	@mkdir -p $(INC_DIR)/GNWA
	+ cd $(ALGO_DIR)/gssw && cp -r src/*.h $(INC_DIR)/gssw && cp -r src/simde $(INC_DIR)/gssw
	+ cd $(ALGO_DIR)/gwfa && cp -r *.h $(INC_DIR)/gwfa
	+ cd $(ALGO_DIR)/ksw2 && cp -r *.h $(INC_DIR)/ksw2
	+ cd $(ALGO_DIR)/csswl/src && cp -r *.h $(INC_DIR)/csswl
	+ cd $(ALGO_DIR)/edlib && cp -r ./edlib/include/edlib.h $(INC_DIR)/edlib
	+ cd $(ALGO_DIR)/abPOA && cp -r ./src/*.h $(INC_DIR)/abPOA
	+ cd $(ALGO_DIR)/vargas && cp -r ./include/* $(INC_DIR)/vargas && cp -r ./cxxopts/src/*.hpp $(INC_DIR)/vargas
	+ cd $(ALGO_DIR)/vargas/htslib && cp -r ./*.h $(INC_DIR)/vargas && cp -r ./htslib/*.h $(INC_DIR)/vargas/htslib
	+ cd $(ALGO_DIR)/vargas/doctest && cp -r ./doctest/*.h $(INC_DIR)/vargas
	+ cd $(ALGO_DIR)/GNWA && cp -r src/*.h $(INC_DIR)/GNWA
	cp $(HEADERS) $(INC_DIR)
	cp $(TEST_HEADERS) $(INC_DIR)
	cp -r $(ALGO_HEADERS) $(INC_DIR)/algorithms

# Linking step for the executable
$(EXE): $(LIBS) $(OBJS) $(MAIN_OBJ) | prepare_headers
	@mkdir -p $(BIN_DIR)
	$(CXX) $(LDFLAGS) -o $@ $(filter %.o, $^) -I$(INC_DIR) $(LIB_FLAGS)

$(TEST_EXE): $(LIBS) $(OBJS) $(TEST_OBJS) | prepare_headers
	@mkdir -p $(BIN_DIR)
	$(CXX) $(LDFLAGS) -o $@ $(filter %.o, $^) -I$(INC_DIR) $(LIB_FLAGS)

# Compile src dir
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | prepare_headers
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $< -I$(INC_DIR)

# Compile main file
$(OBJ_DIR)/main.o: $(SRC_DIR)/main.cpp | prepare_headers
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $< -I$(INC_DIR)

# Compile test dir
$(OBJ_DIR)/%.o: $(TEST_DIR)/%.cpp | prepare_headers
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $< -I$(INC_DIR)

# Compile algorithm dir
$(OBJ_DIR)/%.o: $(ALGO_DIR)/%.cpp | prepare_headers
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $< -I$(INC_DIR)

# Compile deps dir
$(OBJ_DIR)/%.o: $(DEPS_DIR)/



# Create gssw library
$(LIB_DIR)/libgssw.a: prepare_headers
	@mkdir -p $(LIB_DIR)
	+ cd $(ALGO_DIR)/gssw && $(MAKE) && cp -r lib/* $(CWD)/lib

# Create gwfa library
$(LIB_DIR)/libgwfa.a: prepare_headers
	@mkdir -p $(LIB_DIR)
	+ cd $(ALGO_DIR)/gwfa && $(MAKE) && ar rcs libgwfa.a gfa-base.o gfa-io.o gfa-sub.o gwf-ed.o kalloc.o && cp libgwfa.a $(LIB_DIR)

# Create ksw2 library
$(LIB_DIR)/libksw2.a: prepare_headers
	@mkdir -p $(LIB_DIR)
	+ cd $(ALGO_DIR)/ksw2 && $(MAKE) && ar rcs libksw2.a *.o && cp libksw2.a $(LIB_DIR)

# Create csswl library
$(LIB_DIR)/libcsswl.a: prepare_headers
	@mkdir -p $(LIB_DIR)
	+ cd $(ALGO_DIR)/csswl/src && $(MAKE) && ar rcs libcsswl.a *.o && cp libcsswl.a $(LIB_DIR)

# Create edLib library
$(LIB_DIR)/libedlib.a: prepare_headers
	@mkdir -p $(LIB_DIR)
	+ cd $(ALGO_DIR)/edlib && cd build && cmake -D CMAKE_BUILD_TYPE=Release .. && make && cp lib/libedlib.a $(LIB_DIR)

# Create abPOA library
$(LIB_DIR)/libabpoa.a: prepare_headers
	@mkdir -p $(LIB_DIR)
	+ cd $(ALGO_DIR)/abPOA && $(MAKE) && cp -r lib/libabpoa.a $(LIB_DIR)

# Create vargas library
$(LIB_DIR)/libvargas_lib.a: prepare_headers
	@mkdir -p $(LIB_DIR)
	@mkdir -p $(ALGO_DIR)/vargas/build
	+ cd $(ALGO_DIR)/vargas && cd htslib && autoconf && autoheader && $(MAKE) && cd .. && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_AVX512BW_GCC=ON -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc .. && $(MAKE) && cd .. && cp -r lib/libvargas_lib.a $(LIB_DIR)

# Create GNWA library
$(LIB_DIR)/libgnwa.a: prepare_headers
	@mkdir -p $(LIB_DIR)
	@mkdir -p $(ALGO_DIR)/GNWA
	+ cd $(ALGO_DIR)/GNWA && $(MAKE) lib && cp -r lib/libgnwa.a $(LIB_DIR)


run_tests: $(TEST_EXE)
	$(shell $(TEST_EXE))

# Clean up
clean:
	@rm -rf $(OBJ_DIR) $(INC_DIR) $(BIN_DIR) $(LIB_DIR)
	cd $(ALGO_DIR)/gssw && $(MAKE) clean
	cd $(ALGO_DIR)/gwfa && $(MAKE) clean && rm -f libgwfa.a
	cd $(ALGO_DIR)/ksw2 && $(MAKE) clean && rm -f libksw2.a
	cd $(ALGO_DIR)/csswl/src && $(MAKE) clean && rm -f libcsswl.a
	cd $(ALGO_DIR)/edlib && $(MAKE) clean
	cd $(ALGO_DIR)/abPOA && $(MAKE) clean
	cd $(ALGO_DIR)/vargas && cd build && $(MAKE) clean && cd .. && rm -rf bin && rm -rf lib && cd htslib && $(MAKE) clean
	cd $(ALGO_DIR)/GNWA && $(MAKE) clean

.PHONY: all clean prepare_headers