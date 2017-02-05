all: famsa famsa-gpu

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
CFLAGS	= -Wall -O3 -m64 -std=c++11 -pthread -I $(LIBS_DIR) -I $(LIBS_LINUX_DIR) 
CFLAGS_AVX = $(CFLAGS) -mavx ${ABI_FLAG} -mpopcnt -funroll-loops
CFLAGS_AVX2 = $(CFLAGS) -mavx2 ${ABI_FLAG} -mpopcnt -funroll-loops
CLINK	= -lm -O3 -std=c++11 -pthread 

OPENCL_OBJS := 	opencl_utils/hardware/Buffer.o \
	opencl_utils/hardware/DeviceInfo.o \
	opencl_utils/hardware/DeviceWrapper.o \
	opencl_utils/hardware/OpenCl.o \
	opencl_utils/hardware/OwnProgram.o \
	opencl_utils/kernel_repository/KernelFactory.o \
	opencl_utils/kernel_repository/KernelRepository.o \

core/lcsbp_classic.o : core/lcsbp_classic.cpp
	$(CC) $(CFLAGS) -c core/lcsbp_classic.cpp -o $@

core/lcsbp_avx.o : core/lcsbp_avx.cpp
	$(CC) $(CFLAGS_AVX) -c core/lcsbp_avx.cpp -o $@

core/lcsbp_avx2.o : core/lcsbp_avx2.cpp
	$(CC) $(CFLAGS_AVX2) -c core/lcsbp_avx2.cpp -o $@

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

ifeq ($(UNAME_S),Linux)
famsa-gpu: famsa_gpu/famsa_gpu.o \
	core/input_file.o \
	core/msa.o \
	core/output_file.o \
	core/profile.o \
	core/sequence.o \
	core_gpu/gpumsa.o	\
	core/lcsbp_classic.o \
	core/lcsbp_avx.o \
	core/lcsbp_avx2.o \
	core/queues.o \
	core/timer.o \
	core/upgma.o \
	core/NewickTree.o \
	libs/instrset_detect.o \
	opencl_utils/common/Log.o \
	opencl_utils/common/mathex.o \
	opencl_utils/common/StatisticsProvider.o \
	$(OPENCL_OBJS) 
	$(CC) $(CLINK) -o $@ famsa_gpu/famsa_gpu.o \
	core/input_file.o \
	core/msa.o \
	core/output_file.o \
	core/profile.o \
	core/sequence.o \
	core_gpu/gpumsa.o	\
	core/lcsbp_classic.o \
	core/lcsbp_avx.o \
	core/lcsbp_avx2.o \
	core/queues.o \
	core/timer.o \
	core/upgma.o \
	core/NewickTree.o \
	libs/instrset_detect.o \
	opencl_utils/common/Log.o \
	opencl_utils/common/mathex.o \
	opencl_utils/common/StatisticsProvider.o \
	$(OPENCL_OBJS) \
	$(LIBS_LINUX_DIR)/${ASM_LIB} \
	$(LIBS_LINUX_DIR)/CL/libOpenCL.a
endif

famsa: famsa_cpu/famsa_cpu.o \
	core/input_file.o \
	core/msa.o \
	core/output_file.o \
	core/profile.o \
	core/sequence.o \
	core/lcsbp_classic.o \
	core/lcsbp_avx.o \
	core/lcsbp_avx2.o \
	core/queues.o \
	core/timer.o \
	core/upgma.o \
	core/NewickTree.o \
	libs/instrset_detect.o
	$(CC) $(CLINK) -o $@ famsa_cpu/famsa_cpu.o \
	core/input_file.o \
	core/msa.o \
	core/output_file.o \
	core/profile.o \
	core/sequence.o \
	core/lcsbp_classic.o \
	core/lcsbp_avx.o \
	core/lcsbp_avx2.o \
	core/queues.o \
	core/timer.o \
	core/upgma.o \
	core/NewickTree.o \
	libs/instrset_detect.o \
	$(LIBS_LINUX_DIR)/${ASM_LIB}

clean:
	-rm core/*.o
	-rm libs/*.o
	-rm core_gpu/*.o
	-rm opencl_utils/kernel_repository/*.o
	-rm opencl_utils/hardware/*.o
	-rm opencl_utils/common/*.o
	-rm famsa
	-rm famsa-gpu
