all: famsa

## USER'S OPTIONS
STATIC_LINK = false
NO_AVX = false
NO_AVX2 = false
USE_CPP20 = true

####################

LIBS_DIR = libs
LIBS_LINUX_DIR = libs-linux
ABI_FLAG = -fabi-version=0 

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	ABI_FLAG =
endif

ifeq ($(USE_CPP20), true)
	NO_ATOMIC_WAIT = false
	CPP_STD=c++2a
	ATOMIC_FLAGS = 
else
	NO_ATOMIC_WAIT = false
	CPP_STD=c++14
	ATOMIC_FLAGS = -DNO_ATOMIC_WAIT
endif
 
CC 	= g++

ifeq ($(STATIC_LINK), true) 
	CFLAGS	= -Wall -O3 -msse4 -m64 $(ATOMIC_FLAGS) -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=$(CPP_STD) -I $(LIBS_DIR) -I $(LIBS_LINUX_DIR) 
	CLINK	= -lm -static -O3 -msse4 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=$(CPP_STD)
else
	CFLAGS	= -Wall -O3 -msse4 -m64 $(ATOMIC_FLAGS) -std=$(CPP_STD) -pthread -I $(LIBS_DIR) -I $(LIBS_LINUX_DIR)
	CLINK	= -lm -O3 -msse4 -std=$(CPP_STD) -pthread 
endif
 
CFLAGS_AVX = $(CFLAGS) -mavx ${ABI_FLAG} -mpopcnt -funroll-loops
CFLAGS_AVX2 = $(CFLAGS) -mavx2 ${ABI_FLAG} -mpopcnt -funroll-loops


COMMON_OBJS := src/msa.o \
	src/tree/AbstractTreeGenerator.o \
	src/tree/Clustering.o \
	src/tree/GuideTree.o \
	src/tree/NeighborJoining.o \
	src/tree/NewickParser.o \
	src/tree/FastTree.o \
	src/tree/SingleLinkage.o \
	src/tree/UPGMA.o \
	src/utils/pooled_threads.o \
	src/utils/timer.o \
	src/utils/log.o \
	src/core/io_service.o \
	src/core/profile.o \
	src/core/profile_par.o \
	src/core/profile_seq.o \
	src/core/sequence.o \
	src/core/queues.o 
	
src/lcs/lcsbp_classic.o : src/lcs/lcsbp_classic.cpp
	$(CC) $(CFLAGS) -c src/lcs/lcsbp_classic.cpp -o $@

ifeq ($(NO_AVX), true) 
LCS_OBJS := src/lcs/lcsbp.o \
	src/lcs/lcsbp_classic.o
UTILS_OBJS := src/utils/utils.o \
	src/utils/utils.o

src/lcs/lcsbp.o : src/lcs/lcsbp.cpp
	$(CC) $(CFLAGS) -DNO_AVX -c src/lcs/lcsbp.cpp -o $@
src/utils/utils.o : src/utils/utils.cpp
	$(CC) $(CFLAGS) -DNO_AVX -c src/utils/utils.cpp -o $@

else 
ifeq ($(NO_AVX2), true)
LCS_OBJS := src/lcs/lcsbp.o \
	src/lcs/lcsbp_classic.o \
	src/lcs/lcsbp_avx.o \
	src/lcs/lcsbp_avx_intr.o
UTILS_OBJS := src/utils/utils.o \
	src/utils/utils_avx.o 

src/lcs/lcsbp.o : src/lcs/lcsbp.cpp
	$(CC) $(CFLAGS) -DNO_AVX2 -c src/lcs/lcsbp.cpp -o $@
src/lcs/lcsbp_avx_intr.o : src/lcs/lcsbp_avx_intr.cpp
	$(CC) $(CFLAGS_AVX) -c src/lcs/lcsbp_avx_intr.cpp -o $@

src/utils/utils.o : src/utils/utils.cpp
	$(CC) $(CFLAGS_AVX) -c src/utils/utils.cpp -o $@
src/utils/utils_avx.o : src/utils/utils_avx.cpp
	$(CC) $(CFLAGS_AVX) -c src/utils/utils_avx.cpp -o $@
else
LCS_OBJS := src/lcs/lcsbp.o \
	src/lcs/lcsbp_classic.o \
	src/lcs/lcsbp_avx_intr.o \
	src/lcs/lcsbp_avx2_intr.o

UTILS_OBJS := src/utils/utils.o \
	src/utils/utils_avx.o \
	src/utils/utils_avx2.o

src/lcs/lcsbp.o : src/lcs/lcsbp.cpp
	$(CC) $(CFLAGS) -c src/lcs/lcsbp.cpp -o $@
src/lcs/lcsbp_avx_intr.o : src/lcs/lcsbp_avx_intr.cpp
	$(CC) $(CFLAGS_AVX) -c src/lcs/lcsbp_avx_intr.cpp -o $@
src/lcs/lcsbp_avx2_intr.o : src/lcs/lcsbp_avx2_intr.cpp
	$(CC) $(CFLAGS_AVX2) -c src/lcs/lcsbp_avx2_intr.cpp -o $@

src/utils/utils.o : src/utils/utils.cpp
	$(CC) $(CFLAGS) -c src/utils/utils.cpp -o $@
src/utils/utils_avx.o : src/utils/utils_avx.cpp
	$(CC) $(CFLAGS_AVX) -c src/utils/utils_avx.cpp -o $@
src/utils/utils_avx2.o : src/utils/utils_avx2.cpp
	$(CC) $(CFLAGS_AVX2) -c src/utils/utils_avx2.cpp -o $@
endif
endif

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

famsa: src/famsa.o $(COMMON_OBJS) $(LCS_OBJS) $(UTILS_OBJS)
	$(CC) $(CLINK) -o $@ src/famsa.o $(COMMON_OBJS) $(LCS_OBJS) $(UTILS_OBJS)

clean:
	-rm src/core/*.o
	-rm src/lcs/*.o
	-rm src/tree/*.o
	-rm src/utils/*.o
	-rm src/*.o
	-rm famsa

