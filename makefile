all: famsa famsa-gpu

## USER'S OPTIONS
STATIC_LINK = false
NO_AVX = false
NO_AVX2 = false
NO_GPU = false

####################

LIBS_DIR = libs
LIBS_LINUX_DIR = libs-linux
ASM_LIB = libaelf64.a
ABI_FLAG = -fabi-version=0 

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	ASM_LIB = libamac64.a
	ABI_FLAG =
	NO_GPU = true 	
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

OPENCL_OBJS := 	 core_gpu/gpumsa.o \
	opencl_utils/hardware/Buffer.o \
	opencl_utils/hardware/DeviceInfo.o \
	opencl_utils/hardware/DeviceWrapper.o \
	opencl_utils/hardware/OpenCl.o \
	opencl_utils/hardware/OwnProgram.o \
	opencl_utils/kernel_repository/KernelFactory.o \
	opencl_utils/kernel_repository/KernelRepository.o 


COMMON_OBJS := core/input_file.o \
	core/msa.o \
	core/output_file.o \
	core/profile.o \
	core/sequence.o \
	core/queues.o \
	core/timer.o \
	core/upgma.o \
	core/NewickTree.o \
	libs/instrset_detect.o \
	opencl_utils/common/Log.o \
	opencl_utils/common/mathex.o \
	opencl_utils/common/StatisticsProvider.o 

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


ifeq ($(NO_GPU),false)
famsa-gpu: famsa_gpu/famsa_gpu.o $(COMMON_OBJS) $(LCS_OBJS) $(OPENCL_OBJS) 
	$(CC) $(CLINK) -o $@ famsa_gpu/famsa_gpu.o $(COMMON_OBJS) $(LCS_OBJS) $(OPENCL_OBJS) $(LIBS_LINUX_DIR)/${ASM_LIB} -lOpenCL
endif

famsa: famsa_cpu/famsa_cpu.o $(COMMON_OBJS) $(LCS_OBJS)
	$(CC) $(CLINK) -o $@ famsa_cpu/famsa_cpu.o $(COMMON_OBJS) $(LCS_OBJS) $(LIBS_LINUX_DIR)/${ASM_LIB}

clean:
	-rm core/*.o
	-rm libs/*.o
	-rm core_gpu/*.o
	-rm opencl_utils/kernel_repository/*.o
	-rm opencl_utils/hardware/*.o
	-rm opencl_utils/common/*.o
	-rm famsa
	-rm famsa-gpu
