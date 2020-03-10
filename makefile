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


COMMON_OBJS := core/io_service.o \
	core/guide_tree.o \
	core/log.o \
	core/msa.o \
	core/profile.o \
	core/sequence.o \
	core/queues.o \
	core/timer.o \
	core/NewickTree.o \
	libs/instrset_detect.o \
	
core/lcsbp_classic.o : core/lcsbp_classic.cpp
	$(CC) $(CFLAGS) -c core/lcsbp_classic.cpp -o $@

ifeq ($(NO_AVX), true) 
LCS_OBJS := core/lcsbp.o \
	core/lcsbp_classic.o

core/lcsbp.o : core/lcsbp.cpp
	$(CC) $(CFLAGS) -DNO_AVX -c core/lcsbp.cpp -o $@

else 
ifeq ($(NO_AVX2), true)
 
LCS_OBJS := core/lcsbp.o \
	core/lcsbp_classic.o \
	core/lcsbp_avx.o

core/lcsbp.o : core/lcsbp.cpp
	$(CC) $(CFLAGS) -DNO_AVX2 -c core/lcsbp.cpp -o $@
core/lcsbp_avx.o : core/lcsbp_avx.cpp
	$(CC) $(CFLAGS_AVX) -c core/lcsbp_avx.cpp -o $@
else
LCS_OBJS := core/lcsbp.o \
	core/lcsbp_classic.o \
	core/lcsbp_avx.o \
	core/lcsbp_avx2.o

core/lcsbp.o : core/lcsbp.cpp
	$(CC) $(CFLAGS) -c core/lcsbp.cpp -o $@
core/lcsbp_avx.o : core/lcsbp_avx.cpp
	$(CC) $(CFLAGS_AVX) -c core/lcsbp_avx.cpp -o $@
core/lcsbp_avx2.o : core/lcsbp_avx2.cpp
	$(CC) $(CFLAGS_AVX2) -c core/lcsbp_avx2.cpp -o $@
endif
endif

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

famsa: famsa_cpu/famsa_cpu.o $(COMMON_OBJS) $(LCS_OBJS)
	$(CC) $(CLINK) -o $@ famsa_cpu/famsa_cpu.o $(COMMON_OBJS) $(LCS_OBJS) $(LIBS_LINUX_DIR)/${ASM_LIB}

clean:
	-rm core/*.o
	-rm libs/*.o
	-rm famsa_cpu/*.o
	-rm famsa

