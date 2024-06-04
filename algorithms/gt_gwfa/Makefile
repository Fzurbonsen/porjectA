# Compiler
CC = gcc
CFLAGS = -g -Wall -O3 #-fsanitize=address

# Directories
SRC_DIR := src
OBJ_DIR := obj
INC_DIR := include
LIB_DIR := lib
BIN_DIR := bin

# Deps directories
GSSW_DIR := deps/gssw
GWFA_DIR := deps/gwfa
EDLIB_DIR := deps/edlib

# Current working directory
CWD := $(shell pwd)

# Library Flags
LIB_FLAGS = -lz -lgssw -lgwfa -lm -ledlib -lstdc++ #-fsanitize=address
LDFLAGS = -L/$(CWD)/$(LIB_DIR)
# Include Flags
INCLUDE_FLAGS := -I$(CWD)/$(INC_DIR) -isystem $(CWD)/$(INC_DIR) -I$(CWD)/$(SRC_DIR)

# Object files
OBJ := gt_gwfa.o

# Final executable
EXE := test

# Phonies
.PHONY:all bin libs clean
.NOTPARALLEL:all bin libs

all:bin libs
bin:$(BIN_DIR)/$(EXE)
libs:$(LIB_DIR)/libgwfa.a $(LIB_DIR)/libgssw.a $(LIB_DIR)/libgt_gwfa.a



# Rules:
# Rule for binary
$(BIN_DIR)/$(EXE):$(LIB_DIR)/libgwfa.a $(LIB_DIR)/libgssw.a $(LIB_DIR)/libedlib.a $(OBJ_DIR)/$(OBJ) $(SRC_DIR)/main.c
	@mkdir -p $(@D)
	$(CC) $(INCLUDE_FLAGS) $(LDFLAGS) $(CFLAGS) $(SRC_DIR)/main.c $(OBJ_DIR)/$(OBJ) -o $(BIN_DIR)/$(EXE) $(LIB_FLAGS)



# Rule for library libgt_gwfa
$(LIB_DIR)/libgt_gwfa.a: $(OBJ_DIR)/$(OBJ)
	ar rcs $@ $^

# Rule to create objects
$(OBJ_DIR)/$(OBJ):$(SRC_DIR)/gt_gwfa.h $(SRC_DIR)/gt_gwfa.c 
	@mkdir -p $(@D)
	$(CC) $(INCLUDE_FLAGS) $(LDFLAGS) $(CFLAGS) -c -o $@ $(SRC_DIR)/gt_gwfa.c



# Rules to build dependency libraries
$(LIB_DIR)/libgwfa.a: $(wildcard $(GWFA_DIR)/src/*.c) $(wildcard $(GWFA_DIR)/src/*.h)
	@mkdir -p $(@D)
	@mkdir -p $(PWD)/$(INC_DIR)
	@mkdir -p $(PWD)/$(INC_DIR)/gwfa
	+ cd $(GWFA_DIR) && $(MAKE) && cp lib/libgwfa.a $(CWD)/$(LIB_DIR)/ && cp ./*.h $(CWD)/$(INC_DIR)/gwfa/
	
$(LIB_DIR)/libgssw.a: $(GSSW_DIR)/src/gssw.c $(GSSW_DIR)/src/gssw.h
	@mkdir -p $(@D)
	@mkdir -p $(PWD)/$(INC_DIR)
	+ cd $(GSSW_DIR) && $(MAKE) && cp lib/libgssw.a $(CWD)/$(LIB_DIR)/ && cp src/gssw.h $(CWD)/$(INC_DIR)/ && cp -r src/simde $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libedlib.a: $(EDLIB_DIR)/edlib/src/edlib.cpp $(EDLIB_DIR)/edlib/include/edlib.h
	@mkdir -p $(@D)
	@mkdir -p $(PWD)/$(INC_DIR)
	+ cd $(EDLIB_DIR) && cd build && cmake -D CMAKE_BUILD_TYPE=Release .. && make && cp lib/libedlib.a $(CWD)/$(LIB_DIR)/ && cp $(CWD)/$(EDLIB_DIR)/edlib/include/edlib.h $(CWD)/$(INC_DIR)/


# make clean
clean:
	rm -f $(OBJ_DIR)/*.o $(EXE)
	rm -rf $(LIB_DIR) $(INC_DIR) $(BIN_DIR) $(OBJ_DIR)
	cd $(GSSW_DIR) && $(MAKE) clean
	cd $(GWFA_DIR) && $(MAKE) clean
	cd $(EDLIB_DIR) && $(MAKE) clean