ROOT_DIR := /home/elenazaz/prog/BNPlib

PROTO_DIR := $(ROOT_DIR)/mwe/mix/protos
STAN_ROOT_DIR := $(ROOT_DIR)/lib/math
SRC_DIR := $(ROOT_DIR)/mwe/mix/src
SPIKES_DIR := $(SRC_DIR)/spikes
CXX = g++
CFLAGS = \
-fopenmp \
-I$(STAN_ROOT_DIR) \
-I $(STAN_ROOT_DIR)/lib\
-I $(STAN_ROOT_DIR)/lib/eigen_3.3.3/ \
-I$(STAN_ROOT_DIR)/lib/boost_1.72.0/ \
-I$(STAN_ROOT_DIR)/lib/sundials_4.1.0/include \
-I$(STAN_ROOT_DIR)/lib/tbb_2019_U8/include \
-I$(PROTO_DIR)/cpp \
-D_REENTRANT $(shell root-config --cflags)
LDLIBS = \
$(shell pkg-config --libs protobuf) -L$(STAN_ROOT_DIR)/lib/tbb -lpthread -ltbb -Wl,-rpath,"$(STAN_ROOT_DIR)/lib/tbb"
LDFLAGS = -O3 -D_REENTRANT -fopenmp
PROTO_SRCS = $(wildcard $(PROTO_DIR)/cpp/*.cpp)

SPIKES_SRCS = $(wildcard $(SPIKES_DIR)/*.cpp)
OUR_SRCS = $(wildcard $(SRC_DIR)/*.cpp)
SRCS = $(PROTO_SRCS) $(OUR_SRCS)
OBJS = $(subst .cpp,.o, $(SRCS))
SPIKES_EXECS = $(subst .cpp,.out, $(SPIKES_SRCS))
SPIKES_OBJS = $(subst .cpp,.o, $(SPIKES_SRCS))
EXEC = test_main
#info:
	#@echo " Info..." 
	#@echo " ROOT_DIR = $(ROOT_DIR)" 
	#@echo " PROTO_DIR = $(PROTO_DIR)"
	#@echo " SRC_DIR = $(SRC_DIR)"
	#@echo " SPIKES_DIR = $(SPIKES_DIR)"
	#@echo " SOURCES = $(SRCS)"
	#@echo " OBJECTS = $(OBJS)"
	#@echo " EXECS = $(SPIKES_EXECS)"
	#@echo " STAN_ROOT_DIR = $(STAN_ROOT_DIR)"

all: test_main $(SPIKES_EXECS)
#all: $(SPIKES_EXECS) test_main


$(SPIKES_EXECS): %.out: %.o $(OBJS) 
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS) $< $(LDLIBS)
$(SPIKES_OBJS): %.o: %.cpp 
	$(CXX) $(CFLAGS) -c $< -o $@
test_main: test_main.o $(OBJS)
	$(CXX) $(LDFLAGS) -o test_main $(OBJS) test_main.o $(LDLIBS)
test_main.o:
	$(CXX) $(CFLAGS) -c test_main.cpp -o test_main.o
%.o : %.cpp
	$(CXX) $(CFLAGS) -c $< -o $@
-include $(dep)
%.d: %.c 
	$(CXX) $(CFLAGS) $< -MM -MT $(@:.d=.o) >$@


clean:
	rm $(OBJS) $(SPIKES_OBJS) test_main.o
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

