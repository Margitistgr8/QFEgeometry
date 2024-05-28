.SUFFIXES:
.SUFFIXES: .o .cpp
#============================================================
TARGET	=  Tetra

C_SOURCES =  Tetra.cpp 
C_OBJS     =  Tetra.o 
MY_INCLUDES = -I/usr/local/include/eigen 


CCX = g++
CXXFLAGS = -Wall -g -std=c++14 ${MY_INCLUDES}

#============================================================
all: $(TARGET)

.cpp.o:
	$(CCX) -c $(CXXFLAGS) $< -o $@

$(TARGET): $(C_OBJS)
	$(CCX) $(CXXFLAGS) $^ $(LIBDIRS) -o $@


# Implici rules: $@ = target name, $< = first prerequisite name, $^ = name of all prerequisites
#============================================================

NOTES = MathematicaTetra.nb   data/

ALL_SOURCES = Makefile $(C_SOURCES) $(MY_INCLUDES) $(NOTES)

ALL_FILES = $(ALL_SOURCES)  $(NOTES)



clean:
	rm -f $(TARGET) $(C_OBJS) core *.tar *~*~ 

tar: $(ALL_FILES) 
	tar -cvjf $(TARGET).tar $(ALL_FILES)


$(TARGET).ps: $(ALL SOURCES)
	enscript -pcode.ps $(ALL_SOURCES)


