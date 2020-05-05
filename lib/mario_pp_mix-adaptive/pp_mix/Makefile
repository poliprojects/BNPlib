ROOT_DIR := .
STAN_ROOT_DIR := /home/mario/PhD/spatial_lda/lib/math
SRC_DIR := $(ROOT_DIR)/src
SPIKES_DIR := $(SRC_DIR)/spikes
PROTO_DIR := $(ROOT_DIR)/protos

CXX = g++
CFLAGS = \
	-std=c++1y \
	-MMD \
	-I$(STAN_ROOT_DIR) \
	-I$(STAN_ROOT_DIR)/lib/eigen_3.3.3/ \
	-I$(STAN_ROOT_DIR)/lib/boost_1.69.0/  \
	-I$(STAN_ROOT_DIR)/lib/sundials_4.1.0/include  \
	-I$(STAN_ROOT_DIR)/lib/tbb_2019_U8/include  \
	-I$(PROTO_DIR) \
	-D_REENTRANT -fPIC \
	-g # -march=native -msse2 -funroll-loops -ftree-vectorize -fopenmp 
LDLIBS = \
 	$(shell pkg-config --libs protobuf) -L$(STAN_ROOT_DIR)/lib/tbb \
	-lgsl -lgslcblas -lpthread -ltbb -Wl,-rpath,"$(STAN_ROOT_DIR)/lib/tbb"
LDFLAGS = -D_REENTRANT -g # -fopenmp

PROTO_SRCS = $(wildcard $(PROTO_DIR)/cpp/*.cpp)
PG_SRCS = $(wildcard $(PG_DIR)/*.cpp) $(wildcard $(PG_DIR)/include/*.cpp)
SPIKES_SRCS = $(wildcard $(SPIKES_DIR)/*.cpp)
OUR_SRCS = $(wildcard $(SRC_DIR)/*.cpp) \
		   $(wildcard $(SRC_DIR)/precs/*.cpp) \
		   $(wildcard $(SRC_DIR)/jumps/*.cpp) \
		   $(wildcard $(SRC_DIR)/point_process/*.cpp)

PP_SRCS = $(wildcard $(SRC_DIR)/pointprocess/*.cpp)
OUR_SRCS := $(filter-out ./src/python_exports.cpp, $(OUR_SRCS))

SRCS = $(PROTO_SRCS) $(PG_SRCS) $(OUR_SRCS) $(PP_SRCS)
OBJS = $(subst .cpp,.o, $(SRCS))
DEPENDS := $(patsubst %.cpp,%.d,$(SRCS))

SPIKES_EXECS = $(subst .cpp,.out, $(SPIKES_SRCS))
SPIKES_OBJS =  $(subst .cpp,.o, $(SPIKES_SRCS))

info:
	@echo " Info..."
	@echo " ROOT_DIR  = $(ROOT_DIR)"
	@echo " PROTO_DIR = $(PROTO_DIR)"
	@echo " SRC_DIR = $(SRC_DIR)"
	@echo " SPIKES_DIR = $(SPIKES_DIR)"
	@echo " SOURCES = $(SRCS)"
	@echo " OBJECTS = $(OBJS)"
	@echo " EXECS = $(SPIKES_EXECS)"
	@echo " STAN_ROOT_DIR = $(STAN_ROOT_DIR)"

all: generate_pybind $(SPIKES_EXECS)

generate_pybind: $(OBJS)
	$(CXX) -shared $(CFLAGS) `/usr/bin/python3.8 -m pybind11 --includes` \
		src/python_exports.cpp -o \
		pp_mix_cpp`/usr/bin/python3.8-config --extension-suffix` \
		$(OBJS) $(LDLIBS) -fopenmp

$(SPIKES_EXECS): %.out: %.o $(OBJS)
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS) $< $(LDLIBS)

$(SPIKES_OBJS): %.o: %.cpp
		$(CXX) $(CFLAGS) -c $< -o $@

%.o : %.cpp Makefile
	$(CXX) $(CFLAGS) -MMD -MP -c $< -o $@

-include $(DEPENDS)

clean:
	rm $(OBJS) $(SPIKES_OBJS) run_from_file.o

distclean: clean

compile_protos:
	@ mkdir -p $(PROTO_DIR)/cpp;
	@ mkdir -p $(PROTO_DIR)/py;
	@ for filename in $(PROTO_DIR)/*.proto; do \
		protoc --proto_path=$(PROTO_DIR) --python_out=$(PROTO_DIR)/py/ $$filename; \
		protoc --proto_path=$(PROTO_DIR) --cpp_out=$(PROTO_DIR)/cpp/ $$filename; \
	done
	@ for filename in $(PROTO_DIR)/cpp/*.cc; do \
	    mv -- "$$filename" "$${filename%.cc}.cpp"; \
	done

	touch $(PROTO_DIR)/__init__.py
	touch $(PROTO_DIR)/py/__init__.py

	2to3 --output-dir=$(PROTO_DIR)/py/ -W -n $(PROTO_DIR)/py/
