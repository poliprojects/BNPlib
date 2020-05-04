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

LDFLAGS += -O3 -D_REENTRANT -fopenmp
LDLIBS = \
 	$(shell pkg-config --libs protobuf) -L$(STAN_ROOT_DIR)/lib/tbb \
	-lpthread -Wl,-rpath,"$(STAN_ROOT_DIR)/lib/tbb"

EXE = maintest_uni
SRCS_OUTPUT = src/api/output.pb.cc
SRCS = src/api/FileCollector.cpp src/api/MemoryCollector.cpp
OBJS = $(EXE).o $(subst .cc,.o, $(SRCS_OUTPUT)) $(subst .cpp,.o, $(SRCS))

.PHONY: all clean distclean

all: $(EXE) generate_pybind

$(EXE): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(EXE) $(OBJS) $(LDLIBS)

%.h: %.cc %.cpp

-include $(dep)

generate_pybind:
#	$(CXX) -shared $(CFLAGS) `/usr/bin/python3.8 -m pybind11 --includes` \
#		src/python_exports.cpp -o \
#		bnplib`/usr/bin/python3.8-config --extension-suffix` \
#		$(OBJS) $(LDLIBS) -fopenmp

clean:
	$(RM) *.o

distclean: clean
	$(RM) $(EXE)
	$(RM) *~
