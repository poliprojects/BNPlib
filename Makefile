PYTHON_DIR := /usr/bin/python3
STAN_ROOT_DIR := lib/math
EIGEN_DIR := $(STAN_ROOT_DIR)/lib/eigen_3.3.3/
CXX = g++ -fPIC
CXXFLAGS += -fopenmp -D_REENTRANT \
	-I$(STAN_ROOT_DIR) \
	-I$(EIGEN_DIR) \
	-I$(STAN_ROOT_DIR)/lib/boost_1.72.0/ \
	-I$(STAN_ROOT_DIR)/lib/sundials_4.1.0/include \
	-I$(STAN_ROOT_DIR)/lib/tbb_2019_U8/include 
        

LDFLAGS += -O3 -D_REENTRANT -fopenmp
LDLIBS = $(shell pkg-config --libs protobuf) -L$(STAN_ROOT_DIR)/lib/tbb \
	-lpthread -Wl,-rpath,"$(STAN_ROOT_DIR)/lib/tbb"

EXE = maintest_multi
SRCS_OUTPUT = src/api/output.pb.cc
SRCS = src/api/FileCollector.cpp src/api/MemoryCollector.cpp
LIB_OBJS = $(subst .cpp,.o, $(SRCS)) $(subst .cc,.o, $(SRCS_OUTPUT))
OBJS = $(EXE).o $(LIB_OBJS)

.PHONY: all clean distclean

all: $(EXE) pybind_generate

$(EXE): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(EXE) $(OBJS) $(LDLIBS)

%.h: %.cc %.cpp

-include $(dep)

pybind_generate: $(LIB_OBJS)
	$(CXX) -shared $(CXXFLAGS) `$(PYTHON_DIR) -m pybind11 --includes` \
		src/python/exports.cpp -o bnplib`$(PYTHON_DIR)-config \
		--extension-suffix` $(LIB_OBJS) $(LDLIBS)

clean:
	$(RM) *.o
	$(RM) src/api/*.o

distclean: clean
	$(RM) $(EXE)
	$(RM) *~
