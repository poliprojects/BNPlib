STAN_ROOT_DIR := lib/math
CXX = g++
CFLAGS = \
-fopenmp \
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

SRCS = 
OBJS = $(subst .cpp,.o, $(SRCS))
EXEC = main

#info:
	#@echo " Info..." 
	#@echo " ROOT_DIR = $(ROOT_DIR)" 
	#@echo " SRC_DIR = $(SRC_DIR)"
	#@echo " SOURCES = $(SRCS)"
	#@echo " OBJECTS = $(OBJS)"
	#@echo " STAN_ROOT_DIR = $(STAN_ROOT_DIR)"

all: main


main: main.o $(OBJS)
	$(CXX) $(LDFLAGS) -o main $(OBJS) main.o $(LDLIBS)
main.o:
	$(CXX) $(CFLAGS) -c main.cpp -o main.o
%.o : %.cpp
	$(CXX) $(CFLAGS) -c $< -o $@
-include $(dep)
%.d: %.c 
	$(CXX) $(CFLAGS) $< -MM -MT $(@:.d=.o) >$@


clean:
	rm $(OBJS) main.o

distclean: clean
