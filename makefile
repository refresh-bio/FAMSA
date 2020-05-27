all: famsa

## USER'S OPTIONS
STATIC_LINK = false
NO_AVX = false
NO_AVX2 = false

####################

LIBS_DIR = libs
LIBS_LINUX_DIR = libs-linux
ASM_LIB = libaelf64.a
ABI_FLAG = -fabi-version=0 

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	ASM_LIB = libamac64.a
	ABI_FLAG =
endif

 
CC 	= g++

ifeq ($(STATIC_LINK), true) 
	CFLAGS	= -Wall -O3 -m64 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++11 -I $(LIBS_DIR) -I $(LIBS_LINUX_DIR) 
	CLINK	= -lm -static -O3 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++11 
else
	CFLAGS	= -Wall -O3 -m64 -std=c++11 -pthread -I $(LIBS_DIR) -I $(LIBS_LINUX_DIR)
	CLINK	= -lm -O3 -std=c++11 -pthread 
endif
 
CFLAGS_AVX = $(CFLAGS) -mavx ${ABI_FLAG} -mpopcnt -funroll-loops
CFLAGS_AVX2 = $(CFLAGS) -mavx2 ${ABI_FLAG} -mpopcnt -funroll-loops


COMMON_OBJS := src/msa.o \
	src/tree/AbstractTreeGenerator.o \
	src/tree/Clustering.o \
	src/tree/GuideTree.o \
	src/tree/NeighborJoining.o \
	src/tree/NewickParser.o \
	src/tree/PartTree.o \
	src/tree/SingleLinkage.o \
	src/tree/UPGMA.o \
	src/utils/timer.o \
	src/utils/log.o \
	src/core/io_service.o \
	src/core/profile.o \
	src/core/sequence.o \
	src/core/queues.o \
	libs/instrset_detect.o \
	
src/lcs/lcsbp_classic.o : src/lcs/lcsbp_classic.cpp
	$(CC) $(CFLAGS) -c src/lcs/lcsbp_classic.cpp -o $@

ifeq ($(NO_AVX), true) 
LCS_OBJS := src/lcs/lcsbp.o \
	src/lcs/lcsbp_classic.o

src/lcs/lcsbp.o : src/lcs/lcsbp.cpp
	$(CC) $(CFLAGS) -DNO_AVX -c src/lcs/lcsbp.cpp -o $@

else 
ifeq ($(NO_AVX2), true)
 
LCS_OBJS := src/lcs/lcsbp.o \
	src/lcs/lcsbp_classic.o \
	src/lcs/lcsbp_avx.o

src/lcs/lcsbp.o : src/lcs/lcsbp.cpp
	$(CC) $(CFLAGS) -DNO_AVX2 -c src/lcs/lcsbp.cpp -o $@
src/lcs/lcsbp_avx.o : src/lcs/lcsbp_avx.cpp
	$(CC) $(CFLAGS_AVX) -c src/lcs/lcsbp_avx.cpp -o $@
else
LCS_OBJS := src/lcs/lcsbp.o \
	src/lcs/lcsbp_classic.o \
	src/lcs/lcsbp_avx.o \
	src/lcs/lcsbp_avx2.o

src/lcs/lcsbp.o : src/lcs/lcsbp.cpp
	$(CC) $(CFLAGS) -c src/lcs/lcsbp.cpp -o $@
src/lcs/lcsbp_avx.o : src/lcs/lcsbp_avx.cpp
	$(CC) $(CFLAGS_AVX) -c src/lcs/lcsbp_avx.cpp -o $@
src/lcs/lcsbp_avx2.o : src/lcs/lcsbp_avx2.cpp
	$(CC) $(CFLAGS_AVX2) -c src/lcs/lcsbp_avx2.cpp -o $@
endif
endif

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

famsa: src/famsa.o $(COMMON_OBJS) $(LCS_OBJS)
	$(CC) $(CLINK) -o $@ src/famsa.o $(COMMON_OBJS) $(LCS_OBJS) $(LIBS_LINUX_DIR)/${ASM_LIB}

clean:
	-rm src/core/*.o
	-rm src/lcs/*.o
	-rm src/tree/*.o
	-rm src/utils/*.o
	-rm src/*.o
	-rm libs/*.o
	-rm famsa

