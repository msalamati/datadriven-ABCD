# Copyright 2021, Gurobi Optimization, LLC

GRBROOT = /home/bwooding/gurobi912
PLATFORM = linux64
TWOUP    = $(GRBROOT)/$(PLATFORM)/examples
INC      = $(TWOUP)/../include/
CPP      = g++
CARGS    = -m64 -g
CPPLIB   = -L$(TWOUP)/../lib -lgurobi_c++ -lgurobi91
GRBAPP   = DOTNETCore2
DOTNETFRAMEWORK = --framework=netcoreapp2.1
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O -DNDEBUG

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
CUDDLPATH   	= -L$(CUDDPATH)/

%_c++: %_c++.cpp
	$(CPP) $(CARGS) $(CXXFLAGS) -o $@ $< -I$(INC) $(CUDDINC) $(SCOTSINC) $(CPPLIB) -lm $(CUDDLPATH) $(CUDDLIBS)

clean:
	rm -rf *.o *_c *_c++ *.class *.log *.rlp *.lp *.bas *.ilp *.mps *.prm *.dSYM; \
	if [ -d $(GRBAPP) ]; then \
		cd $(GRBAPP); \
		find . ! \( -name dotnetcore2.csproj -o -name . \) -exec rm -rf {} +; \
	fi
