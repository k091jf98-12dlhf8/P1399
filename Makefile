# Makefile for PSSketch project

# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -mavx2 -mbmi -Wall -Wextra -std=c++11 -g #-Wno-unused-parameter -Wno-unused-variable -Wno-implicit-fallthrough -Wno-implicit-fallthrough -Wno-format -Wno-unused-but-set-variable

# Target executable
TARGET = PSSketch

# Source files
SRCS = PSSketch.cpp

# Header files (add other headers if needed)
HEADERS = BOBHash32.h class.h definition.h hash.h LF.h PISketch.h strawman.h para.h bitset.h OO_PE.h CMSketch.h data.h Abstract.h class2.h

# Object files
OBJS = $(SRCS:.cpp=.o)

# Default rule to build the target
all: $(TARGET)

# Rule to build the target executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Rule to build object files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $<

# Rule to clean the build directory
clean:
	rm -f $(OBJS) $(TARGET)

# Phony targets
.PHONY: all clean
