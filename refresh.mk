### REFRESH group macros - v.1.0.16 (2025-03-03)

### Macros for initialization
define INIT_GLOBALS
	$(info *** Initialization of global values ***)
	$(eval INCLUDE_DIRS:=-I.)
	$(eval REFRESH_DIR:=)
	$(eval LIBRARY_FILES:=)
	$(eval LINKER_DIRS:=)
	$(eval C_FLAGS:=)
	$(eval CPP_FLAGS:=)
	$(eval PY_FLAGS:=)
	$(eval DEFINE_FLAGS:=)
	$(eval LINKER_FLAGS:=)
	$(eval CMAKE_OSX_FIX:=)
	$(eval COMPILER_ALLOWED:=)
	$(eval TYPE?=release)
	$(eval PREBUILD_JOBS:=)
	$(eval SRC_DIR:=./src)
	$(eval OBJ_DIR:=./obj)
	$(eval OUT_BIN_DIR:=./bin)
	$(eval AR?=ar)
	$(eval NASM?=nasm)
endef

### Macros for 3rd-party libraries registration
# Add zlib-ng
define ADD_ZLIB_NG
	$(info *** Adding zlib-ng ***)
	$(eval ZLIB_DIR:=$(1))
	$(eval ZLIB_A_DIR:=$(1)/build-g++/zlib-ng)
	$(eval ZLIB_A:=$(ZLIB_A_DIR)/libz.a)
	$(eval INCLUDE_DIRS+=-I$(ZLIB_DIR)/build-g++ -I$(ZLIB_DIR)/build-g++/zlib-ng)
	$(eval LIBRARY_FILES+=$(ZLIB_A))
	$(eval LINKER_DIRS+=-L $(ZLIB_A_DIR))
	$(eval PREBUILD_JOBS+=zlib-ng)

	$(eval zlib-ng: $(ZLIB_A))
	$(eval $(ZLIB_A) : ; \
		cd $(ZLIB_DIR) && cmake $(CMAKE_OSX_FIX) -DCMAKE_CXX_COMPILER=$(CXX) -DCMAKE_C_COMPILER=$(CC) -B build-g++/zlib-ng -S . -DZLIB_COMPAT=ON; cmake --build build-g++/zlib-ng --config Release)
endef

# Propose zlib-ng (to be considered by CHOOSE_...)
define PROPOSE_ZLIB_NG
	$(info *** Proposing zlib-ng ***)
	$(eval ZLIB_DIR:=$(1))
	$(eval ZLIB_A_DIR:=$(1)/build-g++/zlib-ng)
	$(eval ZLIB_A:=$(ZLIB_A_DIR)/libz.a)

	$(eval zlib-ng: $(ZLIB_A))
	$(eval $(ZLIB_A) : ; \
		cd $(ZLIB_DIR) && cmake $(CMAKE_OSX_FIX) -DCMAKE_CXX_COMPILER=$(CXX) -DCMAKE_C_COMPILER=$(CC) -B build-g++/zlib-ng -S . -DZLIB_COMPAT=ON; cmake --build build-g++/zlib-ng --config Release)
endef

# Propose isa-l (to be considered by CHOOSE_...)
define PROPOSE_ISAL
	$(info *** Proposing isal ***)
	$(eval ISAL_DIR:=$(1))
	$(eval ISAL_A_DIR:=$(1)/bin)
	$(eval ISAL_A:=$(1)/bin/isa-l.a)

	$(eval isa-l: $(ISAL_A))
	$(eval $(ISAL_A) : ; \
		cd $(ISAL_DIR) && $(MAKE) -f Makefile.unx)
endef

# Add libdeflate
define ADD_LIBDEFLATE
	$(info *** Adding libdeflate ***)
	$(eval INCLUDE_DIRS+=-I$(1))
	$(eval LIBDEFLATE_DIR:=$(1))
	$(eval LIBDEFLATE_A_DIR:=$(1))
	$(eval LIBDEFLATE_A:=$(1)/build/libdeflate.a)
	$(eval LIBRARY_FILES+=$(LIBDEFLATE_A))
	$(eval LINKER_DIRS+=-L $(LIBDEFLATE_A_DIR))
	$(call TEST_SOFT,cmake)
	$(eval PREBUILD_JOBS+=libdeflate)

	$(eval libdeflate: $(LIBDEFLATE_A))
	$(eval $(LIBDEFLATE_A): ; \
		cd $(LIBDEFLATE_DIR) && cmake $(CMAKE_OSX_FIX) -DCMAKE_CXX_COMPILER=$(CXX) -DCMAKE_C_COMPILER=$(CC) -DLIBDEFLATE_BUILD_SHARED_LIB=OFF -DLIBDEFLATE_BUILD_GZIP=OFF -B build && cmake --build build)
endef

# Add zstd
define ADD_LIBZSTD
	$(info *** Adding libzstd ***)
	$(eval INCLUDE_DIRS+=-I$(1))
	$(eval LIBZSTD_DIR:=$(1))
	$(eval LIBZSTD_A_DIR:=$(1))
	$(eval LIBZSTD_A:=$(1)/lib/libzstd.a)
	$(eval LIBRARY_FILES+=$(LIBZSTD_A))
	$(eval LINKER_DIRS+=-L $(LIBZSTD_A_DIR))
	$(eval PREBUILD_JOBS+=libzstd)

	$(eval libzstd: $(LIBZSTD_A))
	$(eval $(LIBZSTD_A): ; \
		cd $(LIBZSTD_DIR) && $(MAKE))
endef

# Add mimalloc
define ADD_MIMALLOC
	$(info *** Adding mimalloc ***)
	$(eval MIMALLOC_INCLUDE_DIR:=$(1)/include)
	$(eval INCLUDE_DIRS+=-I$(1)/include)
	$(eval MIMALLOC_DIR:=$(1))
	$(eval MIMALLOC_OBJ:=$(1)/mimalloc.o)
	$(eval PREBUILD_JOBS+=mimalloc_obj)

	$(eval mimalloc_obj: $(MIMALLOC_OBJ))
	$(eval $(MIMALLOC_OBJ): ; \
		$(CXX) -DMI_MALLOC_OVERRIDE -O3 -DNDEBUG -fPIC -Wall -Wextra -Wno-unknown-pragmas \
		-fvisibility=hidden -ftls-model=initial-exec -fno-builtin-malloc -c -I $(MIMALLOC_INCLUDE_DIR) \
		$(MIMALLOC_DIR)/src/static.c -o $(MIMALLOC_OBJ))
endef

# Add cdflib
define ADD_CDFLIB
	$(info *** Adding cdflib ***)
	$(eval CDFLIB_INCLUDE_DIR:=$(1))
	$(eval INCLUDE_DIRS+=-I$(1))
	$(eval CDFLIB_DIR:=$(1))
	$(eval CDFLIB_OBJ:=$(1)/cdflib.cpp.o)
	$(eval PREBUILD_JOBS+=cdflib_obj)

	$(eval cdflib_obj: $(CDFLIB_OBJ))
	$(eval $(CDFLIB_OBJ): ; \
		cd $(CDFLIB_DIR) && $(CXX) $(CPP_FLAGS) $(OPTIMIZATION_FLAGS) $(DEFINE_FLAGS) $(INCLUDE_DIRS) -c cdflib.cpp -o cdflib.cpp.o)
endef

# Add REFRESH - parallel queues monitor
define ADD_REFRESH_PARALLEL_QUEUES_MONITOR
	$(info *** Adding refresh - parallel queues monitor ***)
	$(eval REFRESH_PARALLEL_QUEUES_MONITOR_DIR:=$(1)/refresh/parallel_queues/lib/)
	$(eval REFRESH_PARALLEL_QUEUES_MONITOR_OBJ:=$(1)/refresh/parallel_queues/lib/parallel-queues-monitor.cpp.o)
	$(eval PREBUILD_JOBS+=refresh_parallel_queues_monitor_obj)
	$(eval refresh_parallel_queues_monitor_obj: $(REFRESH_PARALLEL_QUEUES_MONITOR_OBJ))
	$(eval $(REFRESH_PARALLEL_QUEUES_MONITOR_OBJ): ; \
		cd $(REFRESH_PARALLEL_QUEUES_MONITOR_DIR) && $(CXX) $(CPP_FLAGS) $(OPTIMIZATION_FLAGS) $(DEFINE_FLAGS) $(INCLUDE_DIRS) -c parallel-queues-monitor.cpp -o parallel-queues-monitor.cpp.o)
endef

# Add RADULS-inplace
define ADD_RADULS_INPLACE
	$(info *** Adding raduls-inplace ***)
	$(eval INCLUDE_DIRS+=-I$(1)/Raduls)
	$(eval RADULS_INPLACE_DIR:=$(1)/Raduls)
	$(eval RADULS_INPLACE_A_DIR:=$(1)/Raduls)
	$(eval RADULS_INPLACE_A:=$(1)/Raduls/libraduls.a)
	$(eval LIBRARY_FILES+=$(RADULS_INPLACE_A))
	$(eval LINKER_DIRS+=-L $(RADULS_INPLACE_A_DIR))
	$(eval PREBUILD_JOBS+=raduls-inplace)

	$(eval raduls-inplace: $(RADULS_INPLACE_A))
	$(eval $(RADULS_INPLACE_A) : ; \
		cd $(RADULS_INPLACE_DIR) && $(MAKE))
endef

# Add igraph
define ADD_IGRAPH
	$(info *** Adding igraph ***)
	$(eval INCLUDE_DIRS+=-I$(1)/include -I$(1)/build/include)
	$(eval IGRAPH_DIR:=$(1))
	$(eval IGRAPH_A_DIR:=$(1)/build/src)
	$(eval IGRAPH_A:=$(IGRAPH_A_DIR)/libigraph.a)
	$(eval LIBRARY_FILES+=$(IGRAPH_A))
	$(eval LINKER_DIRS+=-L $(IGRAPH_A_DIR))
	$(eval IGRAPH_TARGET:=igraph)
	$(call TEST_SOFT,cmake)
	$(call TEST_SOFT,bison)
	$(call TEST_SOFT,flex)
	$(eval PREBUILD_JOBS+=igraph)

	$(eval igraph: $(IGRAPH_A))
	$(eval $(IGRAPH_A): ; \
		$(if $(filter Darwin,$(OS_TYPE)), \
			$(eval IEEE754_DOUBLE_ENDIANNESS_MATCHES_FIX:=-DIEEE754_DOUBLE_ENDIANNESS_MATCHES=TRUE), \
			$(eval IEEE754_DOUBLE_ENDIANNESS_MATCHES_FIX:=) \
		) \
		mkdir -p $(IGRAPH_DIR)/build && cmake $(CMAKE_OSX_FIX) $(IEEE754_DOUBLE_ENDIANNESS_MATCHES_FIX) -DCMAKE_CXX_COMPILER=$(CXX) -DCMAKE_C_COMPILER=$(CC) -S $(IGRAPH_DIR) -B $(IGRAPH_DIR)/build && cmake --build $(IGRAPH_DIR)/build )
endef

# Add SBWT
define ADD_SBWT
	$(info *** Adding SBWT ***)
	$(eval INCLUDE_DIRS+=-I$(1)/include -I$(1)/sdsl-lite/include -I$(1)/SeqIO/include -I$(1)/build/external/sdsl-lite/build/external/libdivsufsort/include/)
	$(eval SBWT_DIR:=$(1))
	$(eval SBWT_A_DIR:=$(1)/build)
	$(eval SBWT_A:=$(SBWT_A_DIR)/libsbwt_static.a)
	$(eval SBWT_SDSL_A:=$(SBWT_A_DIR)/external/sdsl-lite/build/lib/libsdsl.a)
	$(eval SBWT_KMC_CORE_A:=$(SBWT_A_DIR)/external/KMC/build/libkmc_core.a)
	$(eval SBWT_KMC_TOOLS_A:=$(SBWT_A_DIR)/external/KMC/build/libkmc_tools.a)
	$(eval LIBRARY_FILES+=$(SBWT_A) $(SBWT_SDSL_A) $(SBWT_KMC_CORE_A) $(SBWT_KMC_TOOLS_A))
	$(eval LINKER_DIRS+=-L $(SBWT_A_DIR))
	$(eval PREBUILD_JOBS+=sbwt)

	$(eval sbwt: $(SBWT_A) $(SBWT_SDSL_A) $(SBWT_KMC_CORE_A) $(SBWT_KMC_TOOLS_A))
	$(eval $(SBWT_A): ; \
		mkdir -p $(SBWT_DIR)/build && cd $(SBWT_DIR)/build && cmake $(CMAKE_OSX_FIX) -DCMAKE_CXX_COMPILER=$(CXX) -DCMAKE_C_COMPILER=$(CC) .. -DMAX_KMER_LENGTH=32 && $(MAKE) -j)
	$(eval $(SBWT_SDSL_A) : $(SBWT_A))
	$(eval $(SBWT_KMC_CORE_A) : $(SBWT_A))
	$(eval $(SBWT_KMC_TOOLS_A) : $(SBWT_A))
endef

# Add chemfiles
define ADD_CHEMFILES
	$(info *** Adding chemfiles ***)
	$(eval CHEMFILES_DIR:=$(1))
	$(eval CHEMFILES_A_DIR:=$(1)/build)
	$(eval CHEMFILES_A:=$(CHEMFILES_A_DIR)/libchemfiles.a)
	$(eval INCLUDE_DIRS+=-I$(CHEMFILES_DIR)/build/include -I$(CHEMFILES_DIR)/include/)
	$(eval LIBRARY_FILES+=$(CHEMFILES_A))
	$(eval LINKER_DIRS+=-L $(CHEMFILES_A_DIR))
	$(eval PREBUILD_JOBS+=chemfiles)

	$(eval chemfiles: $(CHEMFILES_A))
	$(eval $(CHEMFILES_A) : ; \
		mkdir -p $(CHEMFILES_DIR)/build && cd $(CHEMFILES_DIR)/build && cmake .. && cmake --build . -j)
endef

# Add Pybind11
define ADD_PYBIND11
	$(eval PYBIND11_DIR:=$(1))
	$(eval INCLUDE_DIRS+=-I$(PYBIND11_DIR))
	$(eval INCLUDE_DIRS+=-I $(shell python3 -c "import sysconfig;print(sysconfig.get_paths()['include'])"))
	$(eval PY_EXTENSION_SUFFIX:=$(shell python3-config --extension-suffix))
endef

# Add REFRESH libs
define ADD_REFRESH_LIB
	$(info *** Adding REFRESH libs ***)
	$(eval REFRESH_DIR:=-I$(1))
endef

# Add StatsLib
define ADD_STATS_LIB
	$(info *** Adding StatsLib ***)
	$(eval INCLUDE_DIRS+=-I$(1)/include)
endef

# Add Annoy 
define ADD_ANNOY_LIB
	$(info *** Adding Annoy lib ***)
	$(eval INCLUDE_DIRS+=-I$(1)/include)
endef

# Add hnswlib
define ADD_HNSWLIB
	$(info *** Adding hnswlib ***)
	$(eval INCLUDE_DIRS+=-I$(1))
endef

# Add umappp lib
define ADD_UMAPPP_LIB
	$(info *** Adding UMAPPP lib ***)
	$(eval INCLUDE_DIRS+=-I$(1)/include)
endef

# Add CppIrlba lib
define ADD_CPPIRLBA_LIB
	$(info *** Adding CppIrlba lib ***)
	$(eval INCLUDE_DIRS+=-I$(1)/include)
endef

# Add CppKmeans lib
define ADD_CPPKMEANS_LIB
	$(info *** Adding CppIrlba lib ***)
	$(eval INCLUDE_DIRS+=-I$(1)/include)
endef

# Add aarand lib
define ADD_AARAND_LIB
	$(info *** Adding aarand lib ***)
	$(eval INCLUDE_DIRS+=-I$(1)/include)
endef

# Add knncolle lib
define ADD_KNNCOLLE_LIB
	$(info *** Adding knncolle lib ***)
	$(eval INCLUDE_DIRS+=-I$(1)/include)
endef

# Add Eigen lib
define ADD_EIGEN_LIB
	$(info *** Adding Eigen lib ***)
	$(eval INCLUDE_DIRS+=-I$(1))
endef

# Add NASM
# Idea is that if in Makefile ADD_NASM is called it will override the system one
define ADD_NASM
	$(info *** Adding NASM ***)
	$(eval export PATH := $(abspath $(1)):$(PATH))
	$(if $(filter x86_64,$(ARCH_TYPE)), \
		$(eval dummy_install_nasm:=$(shell \
			if [ ! -f $(1)/nasm ]; then \
				cd $(1) && ./autogen.sh && ./configure && $(MAKE) -j; \
			fi) \
		) \
		$(eval NASM:=$(1)/nasm) \
	)
endef

### Macros configuring compiler/linker flags
# Add os-specific flags for static linking
define SET_STATIC
	$(if $(filter true,$(1)), \
		$(if $(filter Darwin,$(OS_TYPE)), \
			$(eval STATIC_LFLAGS:=-static-libgcc -static-libstdc++ -pthread), \
			$(if $(filter x86_64,$(ARCH_TYPE)), \
				$(eval STATIC_LFLAGS:=-static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive), \
				$(eval STATIC_LFLAGS:=-static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive) \
			)
		)
	)
endef

# Add C, C++ standards
define SET_C_CPP_STANDARDS
	$(eval C_STD:=$(1))
	$(eval CPP_STD:=$(2))
endef

# Define allowed compiler version and type
define SET_COMPILER_VERSION_ALLOWED
	$(eval COMPILER_VERSION_$(strip $(1))_$(strip $(2))_MIN:=$(strip $(3)))
	$(eval COMPILER_VERSION_$(strip $(1))_$(strip $(2))_MAX:=$(strip $(4)))
	$(eval COMPILER_ALLOWED+=COMPILER_VERSION_$(strip $(1))_$(strip $(2)))
endef

# Set source, object and binary directories
define SET_SRC_OBJ_BIN
	$(eval SRC_DIR:=$(1))
	$(eval OBJ_DIR:=$(2))
	$(eval OUT_BIN_DIR:=$(3))
endef

# Set lib directories
define SET_3RD_PARTY
	$(eval 3RD_PARTY_DIR:=$(1))
endef

# *** Utility functions
define LESS_THAN
	$(if $(filter 0,$(shell [ $(1) -lt $(2) ]; echo $$?)),1,0)
endef

define GREATER_THAN
	$(if $(filter 0,$(shell [ $(1) -gt $(2) ]; echo $$?)),1,0)
endef

define IN_RANGE
	$(shell if [ $(COMP) -ge $(MIN) ] && [ $(COMP) -le $(MAX) ]; then echo 1; else echo 0; fi)
endef

define TEST_SOFT
	$(if $(shell command -v $(1) >/dev/null 2>&1 && echo found),, \
		$(error The required software '$(1)' is not installed or not in PATH))
endef

# Check Git commit id and set GIT_COMMIT macro for compilation rule
define SET_GIT_COMMIT
	$(eval GIT_COMMIT:=$(shell git describe --always --dirty))
	$(eval DEFINE_FLAGS:=-DGIT_COMMIT=$(GIT_COMMIT))
endef

# Prepare file variables
define LOAD_FILES
$(eval SRC_$(1)_DIR := $(SRC_DIR)/$(2))                              
$(eval OBJ_$(1)_DIR := $(OBJ_DIR)/$(2))                              
$(eval SRC_$(1) := $(wildcard $(SRC_$(1)_DIR)/*.cpp))         
$(eval OBJ_$(1) := $(patsubst $(SRC_$(1)_DIR)/%.cpp, $(OBJ_$(1)_DIR)/%.cpp.o, $(SRC_$(1)))) 
endef

# Dynamic creation of build rules
define DEFAULT_COMPILE_RULE =
$(OBJ_$(1)_DIR)/%.cpp.o: $(SRC_$(1)_DIR)/%.cpp | prebuild
	@mkdir -p $(OBJ_$(1)_DIR)
	$(CXX) $(CPP_FLAGS) $(OPTIMIZATION_FLAGS) $(DEFINE_FLAGS) $(INCLUDE_DIRS) -MMD -MF $$@.d -c $$< -o $$@
endef

# Dynamic creation of build rules for files in directory
define PREPARE_DEFAULT_COMPILE_RULE
# Source files
$(eval SRC_$(1)_DIR := $(SRC_DIR)/$(2))                              
$(eval OBJ_$(1)_DIR := $(OBJ_DIR)/$(2))                              
$(eval SRC_$(1) := $(wildcard $(SRC_$(1)_DIR)/*.cpp))         
$(eval OBJ_$(1) := $(patsubst $(SRC_$(1)_DIR)/%.cpp, $(OBJ_$(1)_DIR)/%.cpp.o, $(SRC_$(1)))) 
# Compilation rule
$(OBJ_$(1)_DIR)/%.cpp.o: $(SRC_$(1)_DIR)/%.cpp | prebuild
	mkdir -p $(OBJ_$(1)_DIR)
	$(CXX) $(3) $(CPP_FLAGS) $(OPTIMIZATION_FLAGS) $(ARCH_FLAGS) $(DEFINE_FLAGS) $(INCLUDE_DIRS) -MMD -MF $$@.d -c $$< -o $$@
# Dependency files
-include $(OBJ_$(1):.o=.o.d)
endef

# Dynamic creation of build rules for files in directory
define PREPARE_DEFAULT_C_COMPILE_RULE
# Source files
$(eval SRC_$(1)_DIR := $(SRC_DIR)/$(2))
$(eval OBJ_$(1)_DIR := $(OBJ_DIR)/$(2))
$(eval SRC_$(1) := $(wildcard $(SRC_$(1)_DIR)/*.c))
$(eval OBJ_$(1) := $(patsubst $(SRC_$(1)_DIR)/%.c, $(OBJ_$(1)_DIR)/%.c.o, $(SRC_$(1))))
# Compilation rule
$(OBJ_$(1)_DIR)/%.c.o: $(SRC_$(1)_DIR)/%.c | prebuild
	mkdir -p $(OBJ_$(1)_DIR)
	$(CC) $(3) $(C_FLAGS) $(OPTIMIZATION_FLAGS) $(ARCH_FLAGS) $(DEFINE_FLAGS) $(INCLUDE_DIRS) -MMD -MF $$@.d -c $$< -o $$@
# Dependency files
-include $(OBJ_$(1):.o=.o.d)
endef

# Dynamic creation of lib build rules
define DEFAULT_LIB_COMPILE_RULE =
$(OBJ_$(1)_DIR)/%.cpp.o: $(3RD_PARTY_$(1)_DIR)/%.cpp
	@mkdir -p $(OBJ_$(1)_DIR)
	$(CXX) $(CPP_FLAGS) $(OPTIMIZATION_FLAGS) $(DEFINE_FLAGS) $(INCLUDE_DIRS) -MMD -MF $$@.d -c $$< -o $$@
endef

# Dynamic creation of build rules for files in lib directory
define PREPARE_DEFAULT_LIBCPP_COMPILE_RULE
# Source files
$(eval 3RD_PARTY_$(1)_DIR := $(3RD_PARTY_DIR)/$(2))                              
$(eval OBJ_$(1)_DIR := $(OBJ_DIR)/$(2))                              
$(eval 3RD_PARTY_$(1) := $(wildcard $(3RD_PARTY_$(1)_DIR)/*.cpp))         
$(eval OBJ_$(1) := $(patsubst $(3RD_PARTY_$(1)_DIR)/%.cpp, $(OBJ_$(1)_DIR)/%.cpp.o, $(3RD_PARTY_$(1)))) 
# Compilation rule
$(OBJ_$(1)_DIR)/%.cpp.o: $(3RD_PARTY_$(1)_DIR)/%.cpp
	mkdir -p $(OBJ_$(1)_DIR)
	$(CXX) $(3) $(CPP_FLAGS) $(OPTIMIZATION_FLAGS) $(ARCH_FLAGS) $(DEFINE_FLAGS) $(INCLUDE_DIRS) -MMD -MF $$@.d -c $$< -o $$@
# Dependency files
-include $(OBJ_$(1):.o=.o.d)
endef

# Dynamic creation of build rules for files in lib directory
define PREPARE_DEFAULT_LIBC_COMPILE_RULE
# Source files
$(eval 3RD_PARTY_$(1)_DIR := $(3RD_PARTY_DIR)/$(2))                              
$(eval OBJ_$(1)_DIR := $(OBJ_DIR)/$(2))                              
$(eval 3RD_PARTY_$(1) := $(wildcard $(3RD_PARTY_$(1)_DIR)/*.c))
$(eval OBJ_$(1) := $(patsubst $(3RD_PARTY_$(1)_DIR)/%.c, $(OBJ_$(1)_DIR)/%.c.o, $(3RD_PARTY_$(1)))) 
# Compilation rule
$(OBJ_$(1)_DIR)/%.c.o: $(3RD_PARTY_$(1)_DIR)/%.c
	mkdir -p $(OBJ_$(1)_DIR)
	$(CC) $(3) $(CPP_FLAGS) $(OPTIMIZATION_FLAGS) $(ARCH_FLAGS) $(DEFINE_FLAGS) $(INCLUDE_DIRS) -MMD -MF $$@.d -c $$< -o $$@
# Dependency files
-include $(OBJ_$(1):.o=.o.d)
endef

# Check compiler version
define CHECK_COMPILER_VERSION
	$(info *** Checking compiler version ***)
	$(eval COMPILER_DESC:=$(shell command -v $(CXX) >/dev/null 2>&1 && basename $(CXX) | sed 's/-.*//' || echo ""))

	$(if $(COMPILER_DESC),,\
		$(error Compiler does not exist) \
	)

	$(eval COMPILER_VERSION_FULL:=$(shell $(CXX) --version | sed -n '1s/^[^0-9]*\([0-9\.]*\).*$$/\1/p'))
	$(eval COMPILER_VERSION_MAJOR:=$(firstword $(subst ., ,$(COMPILER_VERSION_FULL))))

	$(eval COMPILER_DESC:=$(subst g++,GCC,$(COMPILER_DESC)))
	$(eval COMPILER_DESC:=$(subst clang,CLANG,$(COMPILER_DESC)))

	$(info Compiler: $(COMPILER_DESC))
	$(info Version: $(COMPILER_VERSION_MAJOR))

	$(if $(or $(COMPILER_VERSION_$(COMPILER_DESC)_$(OS_ARCH_TYPE)_MIN),$(COMPILER_VERSION_$(COMPILER_DESC)_$(OS_ARCH_TYPE)_MAX)),\
		,\
		$(error Compiler not supported) \
	)

	$(if $(COMPILER_VERSION_$(COMPILER_DESC)_$(OS_ARCH_TYPE)_MIN), \
		$(if $(filter 1,$(call LESS_THAN,$(COMPILER_VERSION_MAJOR),$(COMPILER_VERSION_$(COMPILER_DESC)_$(OS_ARCH_TYPE)_MIN))), \
			$(error Too low compiler version), \
			$(if $(COMPILER_VERSION_$(COMPILER_DESC)_$(OS_ARCH_TYPE)_MAX), \
				$(if $(filter 1,$(call GREATER_THAN,$(COMPILER_VERSION_MAJOR),$(COMPILER_VERSION_$(COMPILER_DESC)_$(OS_ARCH_TYPE)_MAX))), \
					$(error Too high compiler version) \
				), \
			) 
		), \
	)
endef

# Add type-specifix flags
define SET_FLAGS
	$(if $(filter Linux_x86_64,$(OS_ARCH_TYPE)), \
		$(eval PLATFORM_SPECIFIC_C_FLAGS:=) \
		$(eval PLATFORM_SPECIFIC_CPP_FLAGS:=) \
		$(eval PLATFORM_SPECIFIC_LINKER_FLAGS:=-fabi-version=6), \
		$(if $(filter Linux_aarch64,$(OS_ARCH_TYPE)), \
			$(eval PLATFORM_SPECIFIC_C_FLAGS:=) \
			$(eval PLATFORM_SPECIFIC_CPP_FLAGS:=-ffp-contract=off) \
			$(eval PLATFORM_SPECIFIC_LINKER_FLAGS:=-fabi-version=6), \
			$(if $(filter Darwin_arm64,$(OS_ARCH_TYPE)), \
				$(eval PLATFORM_SPECIFIC_C_FLAGS:=) \
				$(eval PLATFORM_SPECIFIC_CPP_FLAGS:=-ffp-contract=off) \
				$(eval PLATFORM_SPECIFIC_LINKER_FLAGS:=), \
				$(if $(filter Darwin_x86_64,$(OS_ARCH_TYPE)), \
					$(eval PLATFORM_SPECIFIC_C_FLAGS:=) \
					$(eval PLATFORM_SPECIFIC_CPP_FLAGS:=) \
					$(eval PLATFORM_SPECIFIC_LINKER_FLAGS:=) \
				) \
			) \
		) \
	)

	$(eval C_FLAGS+=-std=$(C_STD) -Wall -fPIC -pthread -fpermissive $(PLATFORM_SPECIFIC_C_FLAGS))
	$(eval CPP_FLAGS+=-std=$(CPP_STD) -Wall -fPIC -pthread -fpermissive $(PLATFORM_SPECIFIC_CPP_FLAGS))
	$(eval LINKER_FLAGS+=-lm -lpthread $(PLATFORM_SPECIFIC_LINKER_FLAGS) $(STATIC_LFLAGS))
	$(eval PY_FLAGS:=-Wl,-undefined,dynamic_lookup -shared)


	$(if $(filter release,$(1)), \
		$(eval OPTIMIZATION_FLAGS+=-O3) \
		$(eval C_FLAGS+=) \
		$(eval CPP_FLAGS+= ), \
		$(if $(filter debug,$(1)), \
			$(eval OPTIMIZATION_FLAGS+=-O0 -g) \
			$(eval C_FLAGS+=) \
			$(eval CPP_FLAGS+= ), \
			$(if $(filter ASan,$(1)), \
				$(eval OPTIMIZATION_FLAGS+=-O3 -g) \
				$(eval C_FLAGS+=-fsanitize=address) \
				$(eval CPP_FLAGS+=-fsanitize=address) \
				$(eval LINKER_FLAGS+=-fsanitize=address), \
				$(if $(filter TSan,$(1)), \
					$(eval OPTIMIZATION_FLAGS+=-O3 -g) \
					$(eval C_FLAGS+=-fsanitize=thread) \
					$(eval CPP_FLAGS+=-fsanitize=thread) \
					$(eval LINKER_FLAGS+=-fsanitize=thread -static-libgcc -static-libstdc++), \
					$(if $(filter UBSan,$(1)), \
						$(eval OPTIMIZATION_FLAGS+=-O3 -g) \
						$(eval C_FLAGS+=-fsanitize=undefined) \
						$(eval CPP_FLAGS+=-fsanitize=undefined) \
						$(eval LINKER_FLAGS+=-fsanitize=undefined), \
						$(if $(filter LSan,$(1)), \
							$(eval OPTIMIZATION_FLAGS+=-O3 -g) \
							$(eval C_FLAGS+=-fsanitize=leak) \
							$(eval CPP_FLAGS+=-fsanitize=leak) \
							$(eval LINKER_FLAGS+=-fsanitize=leak), \
							$(if $(filter MSan,$(1)), \
								$(eval OPTIMIZATION_FLAGS+=-O3 -g) \
								$(eval C_FLAGS+=-fsanitize=memory) \
								$(eval CPP_FLAGS+=-fsanitize=memory) \
								$(eval LINKER_FLAGS+=-fsanitize=memory), \
							) \
						) \
					) \
				) \
			) \
		) \
	)

	$(eval CPP_FLAGS_SSE2:=$(CPPFLAGS) -msse2)
	$(eval CPP_FLAGS_SSE4:=$(CPPFLAGS) -msse4)
	$(eval CPP_FLAGS_AVX:=$(CPPFLAGS) -mavx)
	$(eval CPP_FLAGS_AVX2:=$(CPPFLAGS) -mavx2)
	$(eval CPP_FLAGS_AVX512:=$(CPPFLAGS) -mavx512f -mavx512dq)
	$(eval CPP_FLAGS_NEON:=$(CPPFLAGS))

	$(eval INCLUDE_DIRS+=$(REFRESH_DIR))
	$(info Prebuild jobs: $(PREBUILD_JOBS))
prebuild: 
	$(PREBUILD_JOBS)
endef

### Macros checking system and software
# Check for NASM
define CHECK_NASM
	$(eval NASM_VERSION:=$(shell $(NASM) --version 2>/dev/null))
endef

# Choose lib for gzip decompression
define CHOOSE_GZIP_DECOMPRESSION
	$(if $(filter x86_64,$(ARCH_TYPE)), \
    	$(if $(and $(NASM_VERSION),$(ISAL_DIR)), \
			$(eval GZ_TARGET:=isa-l) \
			$(eval PREBUILD_JOBS+=isa-l) \
			$(eval INCLUDE_DIRS+=-I$(ISAL_DIR)/include) ,\
			$(eval GZ_TARGET:=zlib-ng) \
			$(eval PREBUILD_JOBS+=zlib-ng) \
			$(eval INCLUDE_DIRS+=-I$(ZLIB_DIR)/build-g++ -I$(ZLIB_DIR)/build-g++/zlib-ng) \
		), \
		$(eval GZ_TARGET:=zlib-ng) \
		$(eval PREBUILD_JOBS+=zlib-ng) \
		$(eval INCLUDE_DIRS+=-I$(ZLIB_DIR)/build-g++ -I$(ZLIB_DIR)/build-g++/zlib-ng) \
  	)

	$(if $(filter isa-l,$(GZ_TARGET)), \
		$(info ISAL will be used for gzip decompression) \
		$(eval GZ_LIB:=isa-l.a) \
		$(eval LIBRARY_FILES+=$(ISAL_A)) \
		$(eval LINKER_DIRS+=-L $(ISAL_A_DIR))
		$(eval C_FLAGS+=-DREFRESH_USE_IGZIP) \
		$(eval CPP_FLAGS+=-DREFRESH_USE_IGZIP), \
		$(info zlib-ng will be used for gzip decompression) \
		$(eval GZ_LIB:=libz.a) \
		$(eval LIBRARY_FILES+=$(ZLIB_A)) \
		$(eval LINKER_DIRS+=-L $(ZLIB_A_DIR))
		$(eval C_FLAGS+=-DREFRESH_USE_ZLIB) \
		$(eval CPP_FLAGS+=-DREFRESH_USE_ZLIB) \
	)
endef

# Check for OS and architecture
define CHECK_OS_ARCH
	$(if $(MSVC), \
		$(eval OS_TYPE:=windows) \
		$(eval ARCH_TYPE:=x86_64), \
		$(eval OS_TYPE:=$(shell uname -s 2>/dev/null || echo not)) \
		$(eval ARCH_TYPE:=$(shell uname -m 2>/dev/null || echo not)) \
		)

	$(eval OS_ARCH_TYPE:=$(OS_TYPE)_$(ARCH_TYPE))

	$(eval CPU_EXTENSIONS_DEFS:=)

	$(if $(shell grep -q 'avx' /proc/cpuinfo && echo yes),$(eval CPU_EXTENSIONS_DEFS+=-DSIMD_AVX),)
	$(if $(shell grep -q 'avx2' /proc/cpuinfo && echo yes),$(eval CPU_EXTENSIONS_DEFS+=-DSIMD_AVX2),)
	$(if $(shell grep -q 'avx512' /proc/cpuinfo && echo yes),$(eval CPU_EXTENSIONS_DEFS+=-DSIMD_AVX512),)
	$(if $(shell grep -q 'neon' /proc/cpuinfo && echo yes),$(eval CPU_EXTENSIONS_DEFS+=-DSIMD_NEON),)

	$(if $(filter arm8,$(1)), \
		$(eval ARCH_FLAGS:=-march=armv8-a -DARCH_ARM) \
		$(info *** ARMv8 with NEON extensions ***), \
		$(if $(filter m1,$(1)), \
			$(eval ARCH_FLAGS:=-march=armv8.4-a -DARCH_ARM -DSIMD_NEON) \
			$(info *** Apple M1 (or newer) with NEON extensions ***), \
			$(if $(filter sse2,$(1)), \
				$(eval ARCH_FLAGS:=-msse2 -m64 -DARCH_X64 -DSIMD_SSE2) \
				$(info *** x86-64 with SSE2 extensions ***), \
				$(if $(filter avx,$(1)), \
					$(eval ARCH_FLAGS:=-mavx -m64 -DARCH_X64 -DSIMD_AVX) \
					$(info *** x86-64 with AVX extensions ***), \
					$(if $(filter avx2,$(1)), \
						$(eval ARCH_FLAGS:=-mavx2 -m64 -DARCH_X64 -DSIMD_AVX2) \
						$(info *** x86-64 with AVX2 extensions ***), \
						$(if $(filter avx512,$(1)), \
							$(eval ARCH_FLAGS:=-mavx512f -mavx512dq -m64 -DARCH_X64 -DSIMD_AVX512) \
							$(info *** x86-64 with AVX512 extensions ***), \
							$(if $(filter generic,$(1)), \
								$(if $(filter x86_64,$(ARCH_TYPE)), \
									$(eval ARCH_FLAGS:=-DARCH_X64) \
									$(info *** Unspecified platform - using generic compilation for x86_64 ***), \
									$(eval ARCH_FLAGS:=-DARCH_ARM) \
									$(info *** Unspecified platform - using generic compilation for ARM ***)), \
								$(if $(filter x86_64,$(ARCH_TYPE)), \
									$(eval ARCH_FLAGS:=-march=native -DARCH_X64) \
									$(eval ARCH_FLAGS+=$(CPU_EXTENSIONS_DEFS)) \
									$(info *** Unspecified platform - using native compilation for x86_64 ***), \
									$(eval ARCH_FLAGS:=-march=native -DARCH_ARM -DSIMD_NEON) \
									$(eval ARCH_FLAGS+=$(CPU_EXTENSIONS_DEFS)) \
									$(info *** Unspecified platform - using native compilation for ARM ***)))))))))
	
	$(if $(filter Darwin,$(OS_TYPE)), \
		$(eval SDK_PATH:=$(shell $(CXX) -v 2>&1 | grep -- '--with-sysroot' | sed -E 's/.*--with-sysroot=([^ ]+).*/\1/')) \
		$(eval CMAKE_OSX_FIX:=-DCMAKE_OSX_SYSROOT=$(SDK_PATH)) \
	)

	$(if $(filter Darwin,$(OS_TYPE)), \
		$(eval AR_OPT:=-rcs) \
		$(eval PY_AGC_API_CFLAGS:=-Wl,-undefined,dynamic_lookup -fPIC -Wall -shared -std=c++14 -O3), \
		$(eval AR_OPT:=rcs -o) \
		$(eval PY_AGC_API_CFLAGS:=-fPIC -Wall -shared -std=c++14 -O3) \
	)
endef

# Load submodules if necessary
#	$(eval dummy:=$(shell git submodule update --init --recursive))

define INIT_SUBMODULES_FAST
	$(info *** Initialization of submodules (fast) ***)
	$(if $(shell git submodule status | grep '^-'), \
		$(info Initializing and updating submodules...) \
		$(eval dummy:=$(shell git submodule update --init --recursive --jobs=8)), \
		$(info Submodules are already up-to-date or none exist.)
	)
endef

define INIT_SUBMODULES
	$(info *** Initialization of submodules ***)
	$(eval dummy:=$(shell git submodule update --init --recursive --jobs=8))
endef

### Clean library targets
clean-zlib-ng:
	-cd $(ZLIB_DIR) && $(MAKE) -f Makefile.in clean && rm -r build-g++

clean-isa-l:
	-cd $(ISAL_DIR) && $(MAKE) -f Makefile.unx clean

clean-libdeflate:
	-cd $(LIBDEFLATE_DIR) && rm -r build

clean-libzstd:
	-cd $(LIBZSTD_DIR) && $(MAKE) clean

clean-raduls-inplace:
	-cd $(RADULS_INPLACE_DIR) && $(MAKE) clean

clean-igraph:
	-rm -r $(IGRAPH_DIR)/build

clean-mimalloc_obj:
	-rm $(MIMALLOC_OBJ)

clean-cdflib_obj:
	-rm $(CDFLIB_OBJ)

clean-refresh_parallel_queues_monitor_obj:
	-rm $(CDFLIB_OBJ)

clean-sbwt:
	-rm $(SBWT_A)
	-rm $(SBWT_SDSL_A)
	-rm $(SBWT_KMC_CORE_A)
	-rm $(SBWT_KMC_TOOLS_A)
	-rm -r $(SBWT_A_DIR)

### Testing
define show_var
	$(info $(1): $($(1)))
endef

define show_var_opt
	$(if $(1), \
		$(info $(1): $($(1))) \
	)
endef

_testing:
	$(info *** General ***)
	$(call show_var,OS_TYPE)
	$(call show_var,ARCH_TYPE)
	$(call show_var,OS_ARCH_TYPE)
	$(call show_var,ARCH_FLAGS)
	$(call show_var,NASM_VERSION)

	$(info *** Compilers ***)
	$(call show_var,COMPILER_DESC)
	$(call show_var,COMPILER_VERSION_FULL)
	$(call show_var,COMPILER_VERSION_MAJOR)
	$(call show_var,COMPILER_ALLOWED)
	$(foreach desc,\
		$(wordlist 1,$(words $(COMPILER_ALLOWED)),$(COMPILER_ALLOWED)), \
		$(call show_var,$(desc)_MIN) \
		$(call show_var,$(desc)_MAX) \
		)

	$(info *** Main directories ***)
	$(call show_var,INCLUDE_DIRS)
	$(call show_var,LIBRARY_DIRS)

	$(info *** Compiler and linker flags ***)
	$(call show_var,C_STD)
	$(call show_var,CPP_STD)
	$(call show_var,C_FLAGS)
	$(call show_var,CPP_FLAGS)
	$(call show_var,OPTIMIZATION_FLAGS)
	$(call show_var,DEFINE_FLAGS)
	$(call show_var,LINKER_FLAGS)
	$(call show_var,STATIC_LFLAGS)
	$(call show_var,CPP_FLAGS_SSE2)
	$(call show_var,CPP_FLAGS_SSE4)
	$(call show_var,CPP_FLAGS_AVX)
	$(call show_var,CPP_FLAGS_AVX2)
	$(call show_var,CPP_FLAGS_AVX512)
	$(call show_var,CPP_FLAGS_NEON)

	$(info *** Files ***)
	$(call show_var,SRC_DIR)
	$(call show_var,OBJ_DIR)
	$(call show_var,OUT_BIN_DIR)
	$(call show_var,FILES_DEFINED)
	$(foreach item,\
		$(wordlist 1,$(words $(FILES_DEFINED)),$(FILES_DEFINED)), \
		$(call show_var,SRC_$(item)_DIR) \
		$(call show_var,OBJ_$(item)_DIR) \
		$(call show_var,SRC_$(item)) \
		$(call show_var,OBJ_$(item)) \
		)

	$(info *** Libraries ***)
	$(info * gzip decompression *)
	$(call show_var,GZ_TARGET)

	$(info * zlib-ng *)
	$(call show_var_opt,ZLIB_DIR)
	$(call show_var_opt,ZLIB_A_DIR)
	$(call show_var_opt,ZLIB_A)

	$(info * isa-l *)
	$(call show_var_opt,ISAL_DIR)
	$(call show_var_opt,ISAL_A_DIR)
	$(call show_var_opt,ISAL_A)

	$(info * libdeflate *)
	$(call show_var_opt,LIBDEFLATE_DIR)
	$(call show_var_opt,LIBDEFLATE_A_DIR)
	$(call show_var_opt,LIBDEFLATE_A)
	
	$(info * libzstd *)
	$(call show_var_opt,LIBZSTD_DIR)
	$(call show_var_opt,LIBZSTD_A_DIR)
	$(call show_var_opt,LIBZSTD_A)

	$(info * mimalloc *)
	$(call show_var_opt,MIMALLOC_INCLUDE_DIR)
	$(call show_var_opt,MIMALLOC_DIR)
	$(call show_var_opt,MIMALLOC_OBJ)

	$(info * raduls *)
	$(call show_var_opt,RADULS_INPLACE_DIR)
	$(call show_var_opt,RADULS_INPLACE_A_DIR)
	$(call show_var_opt,RADULS_INPLACE_A)

	$(info * igraph *)
	$(call show_var_opt,IGRAPH_DIR)
	$(call show_var_opt,IGRAPH_A_DIR)
	$(call show_var_opt,IGRAPH_A)
	
	$(info * SBWT *)
	$(call show_var_opt,SBWT_DIR)
	$(call show_var_opt,SBWT_A_DIR)
	$(call show_var_opt,SBWT_A)
	$(call show_var_opt,SBWT_SDSL_A)
	$(call show_var_opt,SBWT_KMC_CORE_A)
	$(call show_var_opt,SBWT_KMC_TOOLS_A)

### COMMENTS / EXPLANATIONS
# * -ffp-contract=off - used on ARM architecture, because in SPLASH in some cases we had different results on ARM,
#                       probably exactly the same as here: https://bugs.mysql.com/bug.php?id=82760

##############################
