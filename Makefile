# 1. Compiler and Flags
CXX = g++
# The -I$(INC_DIR) flag tells the compiler to look in the 'include' folder for .h files
CXXFLAGS = -O3 -std=c++17 -Iinclude

# 2. Directories
SRC_DIR = src
INC_DIR = include
OBJ_DIR = build
BIN_DIR = bin

# 3. Object files (Mapped to the build directory)
SHARED_OBJS = $(OBJ_DIR)/Atom_Lookup.o $(OBJ_DIR)/AtomicRadii_Map.o $(OBJ_DIR)/atom.o $(OBJ_DIR)/map.o $(OBJ_DIR)/pdbtovector.o
MAIN_OBJS = $(OBJ_DIR)/main.o $(OBJ_DIR)/internals.o
CLUSTER_OBJS = $(OBJ_DIR)/find_clusters.o $(OBJ_DIR)/cluster.o

# 4. Default target (Builds both executables into the bin/ folder)
all: $(BIN_DIR)/main $(BIN_DIR)/find_clusters

debug: CXXFLAGS = -g -DDEBUG_MODE -std=c++17 -Iinclude
debug: clean all

print: CXXFLAGS = -O3 -DPRINT_MODE -std=c++17 -Iinclude
print: clean all

# 5. Build rules for the executables
$(BIN_DIR)/main: $(MAIN_OBJS) $(SHARED_OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BIN_DIR)/find_clusters: $(CLUSTER_OBJS) $(SHARED_OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

# 6. Generic rule to build .o files from src/%.cpp inside build/
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 7. Rules to create the directories if they don't exist
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# 8. Clean up everything (Deletes the build and bin folders entirely)
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)