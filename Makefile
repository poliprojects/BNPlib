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
OBJS = $(subst .cc,.o, $(SRCS_OUT)) $(subst .cpp,.o, $(SRCS))
PBFLAGS = -pthread -lprotobuf -pthread


.PHONY: all clean distclean

all: main

main: main.o $(OBJS) 
	$(CXX) $(LDFLAGS) -o main $(OBJS) main.o $(LDLIBS) $(PBFLAGS)

main.o: 
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

output.pb.o: output.pb.h

-include $(dep)

clean:
	$(RM) *.o

distclean: clean
	$(RM) $(EXE)
	$(RM) *~
