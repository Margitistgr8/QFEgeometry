CXX=g++
CXXFLAGS=-g -fsanitize=address -std=c++14 -O3 -Wall -Wno-unknown-warning-option \
	-Wno-deprecated-declarations -Wno-sign-compare \
	-I /usr/local/include/eigen -I/usr/local/include/eigen/unsupported

SRC=$(wildcard *.cc)
BIN_DIR=../bin
BIN=$(basename $(SRC))
INC_DIR=../include
INC=$(wildcard $(INC_DIR)/*.h)

.PHONY: clean

all: $(BIN)

%: %.cc $(INC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I $(INC_DIR) $< -o $(BIN_DIR)/$@

clean:
	$(RM) -r $(addprefix $(BIN_DIR)/, $(BIN))
