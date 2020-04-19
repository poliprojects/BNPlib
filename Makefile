STAN_ROOT_DIR := lib/math
CXX = g++
CXXFLAGS += \
-fopenmp \
-I$(STAN_ROOT_DIR) \
-I$(STAN_ROOT_DIR)/lib/eigen_3.3.3/ \
-I$(STAN_ROOT_DIR)/lib/boost_1.72.0/ \
-I$(STAN_ROOT_DIR)/lib/sundials_4.1.0/include \
-I$(STAN_ROOT_DIR)/lib/tbb_2019_U8/include \
-D_REENTRANT

LDFLAGS += -O3 -D_REENTRANT
LDLIBS = \
 	$(shell pkg-config --libs protobuf) -L$(STAN_ROOT_DIR)/lib/tbb \
	-lpthread -Wl,-rpath,"$(STAN_ROOT_DIR)/lib/tbb"


SRCS_OUTPUT = src/api/output.pb.cc
SRCS =
OBJS = main.o $(subst .cc,.o, $(SRCS_OUTPUT)) $(subst .cpp,.o, $(SRCS))
EXE = main

.PHONY: all clean distclean

all: $(EXE)

$(EXE): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(EXE) $(OBJS) $(LDLIBS)


%.h: %.cc

-include $(dep)

clean:
	$(RM) *.o

distclean: clean
	$(RM) main
	$(RM) *~
