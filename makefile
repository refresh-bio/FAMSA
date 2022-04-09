all: deflate famsa 

deflate: 
	$(MAKE) -C libs/libdeflate

####################
LIBS_DIR = libs
LIB_FILES := libs/libdeflate/libdeflate.a

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	GCC_VERSION=8
	ABI_FLAG =
else
	GCC_VERSION = $(shell $(CXX) -dumpversion | cut -f1 -d.)
	ABI_FLAG = -fabi-version=0 
endif


ifeq ($(GCC_VERSION), 5)
$(info *** Detecting g++ version 5 ***)
	CPP_STD=c++11
	DEFINE_FLAGS = -DNO_PROFILE_PAR -DOLD_ATOMIC_FLAG
endif
ifeq ($(GCC_VERSION), 6)
$(info *** Detecting g++ version 6 ***)
	CPP_STD=c++11
	DEFINE_FLAGS = -DNO_PROFILE_PAR -DOLD_ATOMIC_FLAG
endif
ifeq ($(GCC_VERSION), 7)
$(info *** Detecting g++ version 7 ***)
	CPP_STD=c++14
	DEFINE_FLAGS = -DNO_PROFILE_PAR -DOLD_ATOMIC_FLAG
endif
ifeq ($(GCC_VERSION), 8)
$(info *** Detecting g++ version 8 ***)
	CPP_STD=c++2a
	DEFINE_FLAGS = -DOLD_ATOMIC_FLAG
endif
ifeq ($(GCC_VERSION), 9)
$(info *** Detecting g++ version 9 ***)
	CPP_STD=c++2a
	DEFINE_FLAGS =  -DOLD_ATOMIC_FLAG
endif
ifeq ($(GCC_VERSION), 10)
$(info *** Detecting g++ version 10 ***)
	CPP_STD=c++2a
	DEFINE_FLAGS = -DOLD_ATOMIC_FLAG
endif
ifeq ($(GCC_VERSION), 11)
$(info *** Detecting g++ version 11 ***)
	CPP_STD=c++20
	DEFINE_FLAGS = 
endif
ifeq ($(GCC_VERSION), 12)
$(info *** Detecting g++ version 12 ***)
	CPP_STD=c++20
#	DEFINE_FLAGS = -DUSE_NATIVE_BARRIERS
	DEFINE_FLAGS = 
endif
 
 
# Detecting user's options and add flags
ifeq ($(CPU_EXT), none32) 
$(info *** Building 32 bits w/o SIMD extensions ***)
	COMMON_FLAGS :=  
	DEFINE_FLAGS := $(DEFINE_FLAGS) -DNO_AVX
	NO_AVX=1
else ifeq ($(CPU_EXT), none64)
$(info *** Building 64 bits w/o SIMD extensions ***)
	COMMON_FLAGS := -m64
	DEFINE_FLAGS := $(DEFINE_FLAGS) -DNO_AVX
	NO_AVX=1
else ifeq ($(CPU_EXT), sse4)
$(info *** Building with SSE4 extensions***)
	COMMON_FLAGS := -m64 -msse4
	DEFINE_FLAGS := $(DEFINE_FLAGS) -DNO_AVX
	NO_AVX=1
else ifeq ($(CPU_EXT), avx1)
$(info *** Building with AVX extensions***)
	COMMON_FLAGS := -m64 -msse4
	DEFINE_FLAGS := $(DEFINE_FLAGS) -DNO_AVX2
	NO_AVX=2
else
$(info *** Building with AVX and AVX2 extensions***)
	COMMON_FLAGS := -m64 -msse4
endif
 
 

 
 
ifeq ($(STATIC_LINK), true) 
	CFLAGS	= -Wall -Wno-char-subscripts -Wno-attributes -O3 -$(COMMON_FLAGS) $(DEFINE_FLAGS) -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=$(CPP_STD) -I $(LIBS_DIR)
	CLINK	= -lm -static -O3 -msse4 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=$(CPP_STD)
else
	CFLAGS	= -Wall -Wno-char-subscripts -Wno-attributes -O3 $(COMMON_FLAGS) $(DEFINE_FLAGS) -std=$(CPP_STD) -pthread -I $(LIBS_DIR)
	CLINK	= -lm -O3 $(COMMON_FLAGS) -std=$(CPP_STD) -pthread 
endif
 
CFLAGS_AVX = $(CFLAGS) -mavx ${ABI_FLAG} -mpopcnt -funroll-loops
CFLAGS_AVX2 = $(CFLAGS) -mavx2 ${ABI_FLAG} -mpopcnt -funroll-loops



COMMON_OBJS := src/msa.o \
	src/tree/AbstractTreeGenerator.o \
	src/tree/Clustering.o \
	src/tree/DistanceCalculator.o \
	src/tree/FastTree.o \
	src/tree/GuideTree.o \
	src/tree/MSTPrim.o \
	src/tree/NeighborJoining.o \
	src/tree/NewickParser.o \
	src/tree/SingleLinkage.o \
	src/tree/UPGMA.o \
	src/utils/timer.o \
	src/utils/log.o \
	src/core/io_service.o \
	src/core/params.o \
	src/core/profile.o \
	src/core/profile_par.o \
	src/core/profile_seq.o \
	src/core/sequence.o \
	src/core/queues.o 
		
src/lcs/lcsbp_classic.o : src/lcs/lcsbp_classic.cpp
	$(CXX) $(CFLAGS) -c src/lcs/lcsbp_classic.cpp -o $@

ifeq ($(NO_AVX), 1) 
LCS_OBJS := src/lcs/lcsbp.o \
	src/lcs/lcsbp_classic.o
UTILS_OBJS := src/utils/utils.o 

src/lcs/lcsbp.o : src/lcs/lcsbp.cpp
	$(CXX) $(CFLAGS) -c src/lcs/lcsbp.cpp -o $@
src/utils/utils.o : src/utils/utils.cpp
	$(CXX) $(CFLAGS) -c src/utils/utils.cpp -o $@

else ifeq ($(NO_AVX), 2)
LCS_OBJS := src/lcs/lcsbp.o \
	src/lcs/lcsbp_classic.o \
	src/lcs/lcsbp_avx_intr.o
UTILS_OBJS := src/utils/utils.o \
	src/utils/utils_avx.o 

src/lcs/lcsbp.o : src/lcs/lcsbp.cpp
	$(CXX) $(CFLAGS) -c src/lcs/lcsbp.cpp -o $@
src/lcs/lcsbp_avx_intr.o : src/lcs/lcsbp_avx_intr.cpp
	$(CXX) $(CFLAGS_AVX) -c src/lcs/lcsbp_avx_intr.cpp -o $@

src/utils/utils.o : src/utils/utils.cpp
	$(CXX) $(CFLAGS) -c src/utils/utils.cpp -o $@
src/utils/utils_avx.o : src/utils/utils_avx.cpp
	$(CXX) $(CFLAGS_AVX) -c src/utils/utils_avx.cpp -o $@
else
LCS_OBJS := src/lcs/lcsbp.o \
	src/lcs/lcsbp_classic.o \
	src/lcs/lcsbp_avx_intr.o \
	src/lcs/lcsbp_avx2_intr.o

UTILS_OBJS := src/utils/utils.o \
	src/utils/utils_avx.o \
	src/utils/utils_avx2.o

src/lcs/lcsbp.o : src/lcs/lcsbp.cpp
	$(CXX) $(CFLAGS) -c src/lcs/lcsbp.cpp -o $@
src/lcs/lcsbp_avx_intr.o : src/lcs/lcsbp_avx_intr.cpp
	$(CXX) $(CFLAGS_AVX) -c src/lcs/lcsbp_avx_intr.cpp -o $@
src/lcs/lcsbp_avx2_intr.o : src/lcs/lcsbp_avx2_intr.cpp
	$(CXX) $(CFLAGS_AVX2) -c src/lcs/lcsbp_avx2_intr.cpp -o $@

src/utils/utils.o : src/utils/utils.cpp
	$(CXX) $(CFLAGS) -c src/utils/utils.cpp -o $@
src/utils/utils_avx.o : src/utils/utils_avx.cpp
	$(CXX) $(CFLAGS_AVX) -c src/utils/utils_avx.cpp -o $@
src/utils/utils_avx2.o : src/utils/utils_avx2.cpp
	$(CXX) $(CFLAGS_AVX2) -c src/utils/utils_avx2.cpp -o $@
endif


.cpp.o:
	$(CXX) $(CFLAGS) -c $< -o $@

famsa: src/famsa.o $(COMMON_OBJS) $(LCS_OBJS) $(UTILS_OBJS)
	$(CXX) $(CLINK) -o $@ src/famsa.o $(COMMON_OBJS) $(LCS_OBJS) $(UTILS_OBJS) $(LIB_FILES)

clean:
	$(MAKE) clean -C libs/libdeflate
	-rm src/core/*.o
	-rm src/lcs/*.o
	-rm src/tree/*.o
	-rm src/utils/*.o
	-rm src/*.o
	-rm famsa

