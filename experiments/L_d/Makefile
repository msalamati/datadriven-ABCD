#
# compiler
#
CC        = g++
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG

#
# scots
#
SCOTSROOT		= /home/bwooding/scots
SCOTSINC		= -I$(SCOTSROOT)/bdd -I$(SCOTSROOT)/utils

#
# cudd
#
CUDDPATH		=  $(SCOTSROOT)/cudd-3.0.0
CUDDINC 		= -I$(CUDDPATH)
CUDDLIBS		= -lcudd
CUDDLPATH   	= -L$(CUDDPATH)/lib

TARGET = L_d

all: $(TARGET)

%.o:%.cpp
	$(CC) -c $(CXXFLAGS) $(CUDDINC) $(SCOTSINC) $< -o $@

$(TARGET): $(TARGET).o
	$(CC) $(CXXFLAGS) -o $(TARGET) $(TARGET).o $(CUDDLPATH) $(CUDDLIBS)


clean:
	rm  ./$(TARGET)  ./$(TARGET).o
