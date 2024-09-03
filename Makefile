# Compiler and flags
CXX = g++
CXXFLAGS = -g -Wall
#LDFLAGS = -ljson-c  # Link against the json-c library, not used right now, could be useful when parsing json

# Targets and files
TARGET = gen_mutation
SRCS = gen_mutation.cc io.cc utils.cc linkedSequence.cc # Add more source files as needed
OBJS = $(SRCS:.cc=.o)      # Automatically convert .cc files to .o files

# Default target
all: $(TARGET)

# Linking the target
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)
# WITH LDFLAGS: $(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

# Compiling .cc files to .o files
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

# phony commands for clean up and rebuild
.PHONY: clean run rebuild

# Clean up generated files
clean:
	rm -f $(TARGET) $(OBJS)

# Add a separate run target, make run ARGS="arg1 arg2" to pass arguments to the executable., not really used
run: $(TARGET)
	./$(TARGET) $(ARGS)

# Rebuild everything from scratch
rebuild: clean all
