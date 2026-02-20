# 1. Compiler and Flags
CXX = g++
# The -I$(INC_DIR) flag tells the compiler to look in the 'include' folder for .h files
CXXFLAGS = -O3 -std=c++17 -Iinclude

# 2. Directories (MUST BE DEFINED BEFORE TARGETS)
SRC_DIR = src
INC_DIR = include
OBJ_DIR = build
BIN_DIR = bin
RES_DIR = results
TEMP_DIR = temp

# 3. Object files (Mapped to the build directory)
MAIN_OBJS = $(OBJ_DIR)/main.o $(OBJ_DIR)/atom.o $(OBJ_DIR)/Atom_Lookup.o $(OBJ_DIR)/AtomicRadii_Map.o $(OBJ_DIR)/cluster.o $(OBJ_DIR)/internals.o $(OBJ_DIR)/map.o $(OBJ_DIR)/pdbtovector.o $(OBJ_DIR)/pymol.o

# 4. Phony Targets (Commands that are not actual files)
.PHONY: all clean main debug print

# 5. Default and Alias Targets
all: $(BIN_DIR)/allwaters

main: $(BIN_DIR)/allwaters

allwaters: $(BIN_DIR)/allwaters

debug: CXXFLAGS = -g -DDEBUG_MODE -std=c++17 -Iinclude
debug: clean all

print: CXXFLAGS = -O3 -DPRINT_MODE -std=c++17 -Iinclude
print: clean all

# 6. Build rules for the executables
$(BIN_DIR)/allwaters: $(MAIN_OBJS) | $(BIN_DIR) $(RES_DIR) $(TEMP_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

# 7. Generic rule to build .o files from src/%.cpp inside build/
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 8. Rules to create the directories if they don't exist
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(RES_DIR):
	mkdir -p $(RES_DIR)

$(TEMP_DIR):
	mkdir -p $(TEMP_DIR)

# 9. Clean up everything (Deletes the build and bin folders entirely)
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR) $(TEMP_DIR)