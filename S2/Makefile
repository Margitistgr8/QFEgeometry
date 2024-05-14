# Compiler to use
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -g -std=c++14
INCLUDES = -I/usr/local/include -I/usr/local/include/eigen -I/usr/local/include/eigen/unsupported
INCLUDES += -L/usr/local/lib/
# Target executable name
TARGET = affineFlow

# Default target
all: $(TARGET)

# Rule for building the final target
$(TARGET): affineFlow.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(TARGET) affineFlow.o

# Rule for building object files
affineFlow.o: affineFlow.cpp util.h triangle.h 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c affineFlow.cpp

# Clean target
clean:
	rm -f $(TARGET) *.o

