CXX=g++-14

CXXFLAGS=-g -std=c++14 -O3 -fopenmp -Wall -Wno-unknown-warning-option \
	-Wno-deprecated-declarations -Wno-sign-compare \
	-I /usr/local/include/eigen -I/usr/local/include/eigen/unsupported \
	-I/opt/homebrew/Cellar/boost/1.88.0/include \

SRC=$(wildcard *.cc)
BIN_DIR=../bin
BIN=$(basename $(SRC))
INC_DIR=../include
INC=$(wildcard $(INC_DIR)/*.h)


GRP_DIR=$(shell pwd)/../grp
CXXFLAGS+=-DGRP_DIR="\"$(GRP_DIR)\""


.PHONY: clean

all: $(BIN)

%: %.cc $(INC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I $(INC_DIR) $< -o $(BIN_DIR)/$@

clean:
	$(RM) -r $(addprefix $(BIN_DIR)/, $(BIN))
