STAN_ROOT_DIR := lib/math
CXX = g++
CFLAGS = \
-fopenmp \
-I armadillo/include \
-I eigen/eigen-3.3.7/include \
-I$(STAN_ROOT_DIR) \
-I$(STAN_ROOT_DIR)/lib \
-I$(STAN_ROOT_DIR)/lib/eigen_3.3.3/ \
-I$(STAN_ROOT_DIR)/lib/boost_1.72.0/ \
-I$(STAN_ROOT_DIR)/lib/sundials_4.1.0/include \
-I$(STAN_ROOT_DIR)/lib/tbb_2019_U8/include \
-D_REENTRANT $(shell root-config --cflags)
LDLIBS = \
-L$(STAN_ROOT_DIR)/lib/tbb -lpthread -ltbb -Wl,-rpath,"$(STAN_ROOT_DIR)/lib/tbb"
LDFLAGS = -O3 -D_REENTRANT -fopenmp

SRCS = main.cpp Neal8_NNIG.cpp NNIGHierarchy.cpp
OBJS = $(subst .cpp,.o, $(SRCS))
EXEC = main

.PHONY: all clean distclean

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(LDFLAGS) -o main $(OBJS) main.o $(LDLIBS)

main.o: headers.hpp includes.hpp NNIGHierarchy.hpp HypersFixed.hpp \
		Neal8_NNIG.hpp SimpleMixture.hpp

NNIGHierarchy.o: NNIGHierarchy.hpp
HypersFixed.o: HypersFixed.hpp
Neal8_NNIG.o: Neal8_NNIG.hpp
SimpleMixture.o: SimpleMixture.hpp



#info:
	#@echo " Info..." 
	#@echo " ROOT_DIR = $(ROOT_DIR)" 
	#@echo " SRC_DIR = $(SRC_DIR)"
	#@echo " SOURCES = $(SRCS)"
	#@echo " OBJECTS = $(OBJS)"
	#@echo " STAN_ROOT_DIR = $(STAN_ROOT_DIR)"


clean:
	$(RM) *.o

distclean: clean
	$(RM) $(EXE)
	$(RM) *~
