STAN_ROOT_DIR := lib/math
CXX = g++
CXXFLAGS += \
-fopenmp \
-I$(STAN_ROOT_DIR) \
-I$(STAN_ROOT_DIR)/lib \
-I$(STAN_ROOT_DIR)/lib/eigen_3.3.3/ \
-I$(STAN_ROOT_DIR)/lib/boost_1.72.0/ \
-I$(STAN_ROOT_DIR)/lib/sundials_4.1.0/include \
-I$(STAN_ROOT_DIR)/lib/tbb_2019_U8/include \
-D_REENTRANT
LDFLAGS = -O3 -D_REENTRANT -fopenmp
SRCS = NNIGHierarchy.cpp Neal8_NNIG.cpp
OBJS = $(subst .cpp,.o, $(SRCS))

.PHONY: all clean distclean

all: main

main: main.o $(OBJS)
	$(CXX) $(LDFLAGS) -o main $(OBJS) main.o $(LDLIBS)

main.o: includes_main.hpp includes_universal.hpp NNIGHierarchy.hpp \
HypersFixed.hpp Neal8_NNIG.hpp SimpleMixture.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

NNIGHierarchy.o: NNIGHierarchy.hpp HypersFixed.hpp
HypersFixed.o: HypersFixed.hpp
Neal8_NNIG.o: Neal8_NNIG.hpp
SimpleMixture.o: SimpleMixture.hpp

%.h: %.cpp

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

-include $(dep)

clean:
	$(RM) *.o

distclean: clean
	$(RM) $(EXE)
	$(RM) *~
