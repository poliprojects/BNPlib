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
-Ilib/protobuf/src/ \
-Ilib/protobuf/src \
-D_REENTRANT
LDFLAGS += -O3 -D_REENTRANT -fopenmp \
-Llib/protobuf/src/.libs -lprotobuf -pthread

SRCS_OUTPUT = #output.pb.cc
SRCS =
OBJS = main.o $(subst .cc,.o, $(SRCS_OUTPUT)) $(subst .cpp,.o, $(SRCS))

.PHONY: all clean distclean

all: main

main: $(OBJS)
	$(CXX) $(LDFLAGS) -o main $(OBJS)

#output.pb.o: output.pb.h

%.h: %.cc

-include $(dep)

clean:
	$(RM) *.o

distclean: clean
	$(RM) $(EXE)
	$(RM) *~
