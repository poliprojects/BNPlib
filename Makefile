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
-Ilib/protocol_buffers/protobuf-3.11.3/src/ \
-Ilib/protocol_buffers/protobuf-3.11.3/src \
-D_REENTRANT
LDFLAGS = -O3 -D_REENTRANT -fopenmp
SRCS_OUT = output.pb.cc
SRCS =
OBJS = main.o $(subst .cc,.o, $(SRCS_OUT)) $(subst .cpp,.o, $(SRCS))
PBFLAGS = -pthread

.PHONY: all clean distclean

all: main

main: $(OBJS)
	$(CXX) $(LDFLAGS) -o main main.o $(LDLIBS) $(PBFLAGS)

#$(OBJS): output.pb.h HypersFixed.hpp Neal8_NNIG.hpp NNIGHierarchy_imp.hpp \
#Neal8_NNIG_imp.hpp SimpleMixture.hpp NNIGHierarchy.hpp

output.pb.o: output.pb.h

%.h: %.cc

-include $(dep)

clean:
	$(RM) *.o

distclean: clean
	$(RM) $(EXE)
	$(RM) *~
