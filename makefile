all: deflate famsa 

deflate: 
	$(MAKE) -C libs/libdeflate

####################
LIBS_DIR = libs
LIB_FILES := libs/libdeflate/libdeflate.a

UNAME_S := $(shell uname -s)
GCC_VERSION= $(shell $(CXX) -dumpversion | cut -f1 -d.)

ifeq ($(UNAME_S),Darwin)
	ABI_FLAG =
	CLINK_FLAGS =
else
	ABI_FLAG = -fabi-version=0 
	CLINK_FLAGS = -lrt
endif


ifeq ($(GCC_VERSION), 5)
$(info *** Detecting g++ version 5 ***)
	CPP_STD=c++11
	DEFINE_FLAGS = -DNO_PROFILE_PAR -DOLD_ATOMIC_FLAG
else ifeq ($(GCC_VERSION), 6)
$(info *** Detecting g++ version 6 ***)
	CPP_STD=c++11
	DEFINE_FLAGS = -DNO_PROFILE_PAR -DOLD_ATOMIC_FLAG
else ifeq ($(GCC_VERSION), 7)
$(info *** Detecting g++ version 7 ***)
	CPP_STD=c++14
	DEFINE_FLAGS = -DNO_PROFILE_PAR -DOLD_ATOMIC_FLAG
else ifeq ($(GCC_VERSION), 8)
$(info *** Detecting g++ version 8 ***)
	CPP_STD=c++2a
	DEFINE_FLAGS = -DOLD_ATOMIC_FLAG
else ifeq ($(GCC_VERSION), 9)
$(info *** Detecting g++ version 9 ***)
	CPP_STD=c++2a
	DEFINE_FLAGS =  -DOLD_ATOMIC_FLAG
else ifeq ($(GCC_VERSION), 10)
$(info *** Detecting g++ version 10 ***)
	CPP_STD=c++2a
	DEFINE_FLAGS = -DOLD_ATOMIC_FLAG
else ifeq ($(GCC_VERSION), 11)
$(info *** Detecting g++ version 11 ***)
	ifeq ($(UNAME_S),Darwin)
		CPP_STD=c++2a
		DEFINE_FLAGS = -DOLD_ATOMIC_FLAG
	else
		CPP_STD=c++20
		DEFINE_FLAGS = 
	endif
else
$(info *** Detecting g++ version 12 or higher ***)
	ifeq ($(UNAME_S),Darwin)
		CPP_STD=c++2a
		DEFINE_FLAGS = -DOLD_ATOMIC_FLAG
	else
		CPP_STD=c++20
		DEFINE_FLAGS = 
	endif
#	DEFINE_FLAGS = -DUSE_NATIVE_BARRIERS
endif

SIMD_NONE=0
SIMD_AVX1=1
SIMD_AVX2=2
SIMD_AVX512=3
SIMD_NEON=4

 
# Detecting user's options and add flags
ifeq ($(PLATFORM), none) 
$(info *** Unspecified platform w/o extensions ***)
	COMMON_FLAGS :=  
	DEFINE_FLAGS := $(DEFINE_FLAGS) -DSIMD=$(SIMD_NONE)
	SIMD=NONE
else ifeq ($(PLATFORM), arm8)
$(info *** ARMv8 with NEON extensions ***)
	COMMON_FLAGS := -march=armv8-a
	DEFINE_FLAGS := $(DEFINE_FLAGS) -DSIMD=$(SIMD_NEON)
	SIMD=NEON
else ifeq ($(PLATFORM), m1)
$(info *** Apple M1 with NEON extensions ***)
	COMMON_FLAGS := -march=armv8.4-a
	DEFINE_FLAGS := $(DEFINE_FLAGS) -DSIMD=$(SIMD_NEON)
	SIMD=NEON
else ifeq ($(PLATFORM), sse4)
$(info *** x86-64 with SSE4 extensions***)
	COMMON_FLAGS := -msse4
	DEFINE_FLAGS := $(DEFINE_FLAGS) -DSIMD=$(SIMD_NONE)
	SIMD=NONE
else ifeq ($(PLATFORM), avx)
$(info *** x86-64 with AVX extensions***)
	COMMON_FLAGS := -msse4
	DEFINE_FLAGS := $(DEFINE_FLAGS) -DSIMD=$(SIMD_AVX1)
	SIMD=AVX1
else ifeq ($(PLATFORM), native)
$(info *** x86-64 with AVX2 extensions and native architecture ***)
	COMMON_FLAGS := -mavx2 -march=native
	DEFINE_FLAGS := $(DEFINE_FLAGS) -DSIMD=$(SIMD_AVX2)
	SIMD=AVX2
else
$(info *** x86-64 with AVX2 extensions***)
	COMMON_FLAGS := -msse4
	DEFINE_FLAGS := $(DEFINE_FLAGS) -DSIMD=$(SIMD_AVX2)
	SIMD=AVX2
endif
 

# get commit hash
GIT_COMMIT = $(shell git describe --always --dirty) 
DEFINE_FLAGS := $(DEFINE_FLAGS) -DGIT_COMMIT=$(GIT_COMMIT)
 
 
ifeq ($(STATIC_LINK), true) 
	CXXFLAGS	= -Wall -Wno-char-subscripts -Wno-attributes -O3 -$(COMMON_FLAGS) $(DEFINE_FLAGS) -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=$(CPP_STD) -I $(LIBS_DIR)
	CLINK	= -lm -static -O3 -msse4 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=$(CPP_STD)
else
	CXXFLAGS	= -Wall -Wno-char-subscripts -Wno-attributes -O3 $(COMMON_FLAGS) $(DEFINE_FLAGS) -std=$(CPP_STD) -pthread -I $(LIBS_DIR)
	CLINK	= -lm $(CLINK_FLAGS) -O3 $(COMMON_FLAGS) -std=$(CPP_STD) -pthread 
endif
 
CXXFLAGS_AVX = $(CXXFLAGS) -mavx ${ABI_FLAG} -mpopcnt -funroll-loops
CXXFLAGS_AVX2 = $(CXXFLAGS) -mavx2 ${ABI_FLAG} -mpopcnt -funroll-loops
CXXFLAGS_NEON = $(CXXFLAGS) ${ABI_FLAG} -funroll-loops



COMMON_OBJS := src/msa.o \
	src/msa_refinement.o \
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
	src/core/queues.o \
	libs/mimalloc/static.o
		
src/lcs/lcsbp_classic.o : src/lcs/lcsbp_classic.cpp
	$(CXX) $(CXXFLAGS) -c src/lcs/lcsbp_classic.cpp -o $@

ifeq ($(SIMD), NONE) 
LCS_OBJS := src/lcs/lcsbp.o \
	src/lcs/lcsbp_classic.o
UTILS_OBJS := src/utils/utils.o 

src/lcs/lcsbp.o : src/lcs/lcsbp.cpp
	$(CXX) $(CXXFLAGS) -c src/lcs/lcsbp.cpp -o $@
src/utils/utils.o : src/utils/utils.cpp
	$(CXX) $(CXXFLAGS) -c src/utils/utils.cpp -o $@

else ifeq ($(SIMD), AVX1)
LCS_OBJS := src/lcs/lcsbp.o \
	src/lcs/lcsbp_classic.o \
	src/lcs/lcsbp_avx_intr.o
UTILS_OBJS := src/utils/utils.o \
	src/utils/utils_avx.o 

src/lcs/lcsbp.o : src/lcs/lcsbp.cpp
	$(CXX) $(CXXFLAGS) -c src/lcs/lcsbp.cpp -o $@
src/lcs/lcsbp_avx_intr.o : src/lcs/lcsbp_avx_intr.cpp
	$(CXX) $(CXXFLAGS_AVX) -c src/lcs/lcsbp_avx_intr.cpp -o $@

src/utils/utils.o : src/utils/utils.cpp
	$(CXX) $(CXXFLAGS) -c src/utils/utils.cpp -o $@
src/utils/utils_avx.o : src/utils/utils_avx.cpp
	$(CXX) $(CXXFLAGS_AVX) -c src/utils/utils_avx.cpp -o $@
else ifeq ($(SIMD), NEON)
LCS_OBJS := src/lcs/lcsbp.o \
	src/lcs/lcsbp_classic.o \
	src/lcs/lcsbp_neon_intr.o
UTILS_OBJS := src/utils/utils.o \
	src/utils/utils_neon.o 

src/lcs/lcsbp.o : src/lcs/lcsbp.cpp
	$(CXX) $(CXXFLAGS) -c src/lcs/lcsbp.cpp -o $@
src/lcs/lcsbp_neon_intr.o : src/lcs/lcsbp_neon_intr.cpp
	$(CXX) $(CXXFLAGS_NEON) -c src/lcs/lcsbp_neon_intr.cpp -o $@

src/utils/utils.o : src/utils/utils.cpp
	$(CXX) $(CXXFLAGS) -c src/utils/utils.cpp -o $@
src/utils/utils_neon.o : src/utils/utils_neon.cpp
	$(CXX) $(CXXFLAGS_NEON) -c src/utils/utils_neon.cpp -o $@
else
LCS_OBJS := src/lcs/lcsbp.o \
	src/lcs/lcsbp_classic.o \
	src/lcs/lcsbp_avx_intr.o \
	src/lcs/lcsbp_avx2_intr.o

UTILS_OBJS := src/utils/utils.o \
	src/utils/utils_avx.o \
	src/utils/utils_avx2.o

src/lcs/lcsbp.o : src/lcs/lcsbp.cpp
	$(CXX) $(CXXFLAGS) -c src/lcs/lcsbp.cpp -o $@
src/lcs/lcsbp_avx_intr.o : src/lcs/lcsbp_avx_intr.cpp
	$(CXX) $(CXXFLAGS_AVX) -c src/lcs/lcsbp_avx_intr.cpp -o $@
src/lcs/lcsbp_avx2_intr.o : src/lcs/lcsbp_avx2_intr.cpp
	$(CXX) $(CXXFLAGS_AVX2) -c src/lcs/lcsbp_avx2_intr.cpp -o $@

src/utils/utils.o : src/utils/utils.cpp
	$(CXX) $(CXXFLAGS) -c src/utils/utils.cpp -o $@
src/utils/utils_avx.o : src/utils/utils_avx.cpp
	$(CXX) $(CXXFLAGS_AVX) -c src/utils/utils_avx.cpp -o $@
src/utils/utils_avx2.o : src/utils/utils_avx2.cpp
	$(CXX) $(CXXFLAGS_AVX2) -c src/utils/utils_avx2.cpp -o $@
endif


.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

famsa: deflate src/famsa.o $(COMMON_OBJS) $(LCS_OBJS) $(UTILS_OBJS)
	$(CXX) $(CLINK) -o $@ src/famsa.o $(COMMON_OBJS) $(LCS_OBJS) $(UTILS_OBJS) $(LIB_FILES)

clean:
	$(MAKE) clean -C libs/libdeflate
	-rm src/core/*.o
	-rm src/lcs/*.o
	-rm src/tree/*.o
	-rm src/utils/*.o
	-rm src/*.o
	-rm libs/mimalloc/*.o
	-rm famsa

