# gpusieve and its tools (see README.txt).
#
#   make            gpusieve, absort, absplit, ab2keys and gsc2ab
#   make gpusieve   only the GPU siever, which needs the CUDA Toolkit (nvcc)
#   make tools      only absort, absplit, ab2keys and gsc2ab, which need no CUDA
#   make clean
#
# Settings, on the command line (make CUDA_ARCH="70 80" NVCC_CCBIN=g++-13):
#   CUDA_ARCH    compute capabilities to build gpusieve for: "native" (the GPUs of this machine, the default; without a
#                visible GPU nvcc warns and builds for its oldest default instead), "all" (every one this nvcc supports,
#                for a binary to run on other machines), "all-major" (only the X.0 ones, which also
#                run on the X.y GPUs) or a list such as "75 80 86 89 90 120", which also adds PTX for the last one, so that
#                newer GPUs can run it.
#   NVCC         the CUDA compiler (default nvcc)
#   NVCC_THREADS how many compute capabilities nvcc compiles at once (default 0: as many as there are CPUs)
#   NVCC_CCBIN   the host compiler for nvcc, if the default one is newer than the CUDA version supports (CUDA 12.4: GCC 13 at
#                most), for example NVCC_CCBIN=g++-13
#   OPT          optimization flags of the host code (default -O3 -march=native: the programs run on the machine they are
#                built on)
#   CUDA_PTX     the compute capability to embed as PTX with a CUDA_ARCH list (default: the last one in the list), so that
#                GPUs newer than the list can JIT it
#   HAVE_BOINC=1 build gpusieve as a BOINC application: BOINC_DIR is a built BOINC source tree (api/libboinc_api.a and
#                lib/libboinc.a, from ./_autosetup && ./configure --disable-server --disable-client --disable-manager &&
#                make). The host code is then built for any x86-64 (OPT=-O2) with libstdc++ and libgcc linked statically.
#
#   BOINC release build (Rocky 8 image, glibc 2.28):
#     make gpusieve HAVE_BOINC=1 BOINC_DIR=/path/to/boinc CUDA_ARCH="50 52 60 61 70 75 80 86 89 90 100 120" CUDA_PTX=80

CC = gcc
CXX = g++
NVCC ?= nvcc
NVCC_CCBIN ?=
NVCC_THREADS ?= 0
CUDA_ARCH ?= native
HAVE_BOINC ?= 0
BOINC_DIR ?=
ifeq ($(HAVE_BOINC),1)
OPT ?= -O2
else
OPT ?= -O3 -march=native
endif

ifneq ($(filter $(CUDA_ARCH),native all all-major),)
NVCC_ARCH = -arch=$(CUDA_ARCH)
else
CUDA_PTX ?= $(lastword $(CUDA_ARCH))
NVCC_ARCH = $(foreach a,$(CUDA_ARCH),-gencode arch=compute_$(a),code=sm_$(a)) \
            -gencode arch=compute_$(CUDA_PTX),code=compute_$(CUDA_PTX)
endif
NVCCFLAGS = -std=c++17 -O3 $(NVCC_ARCH) -t $(NVCC_THREADS) $(if $(NVCC_CCBIN),-ccbin $(NVCC_CCBIN)) -Xcompiler "$(OPT)" \
            -Wno-deprecated-gpu-targets

ifeq ($(HAVE_BOINC),1)
GLUE_FLAGS = -DHAVE_BOINC -I$(BOINC_DIR) -I$(BOINC_DIR)/api -I$(BOINC_DIR)/lib
GS_LIBS = $(BOINC_DIR)/api/libboinc_api.a $(BOINC_DIR)/lib/libboinc.a -lpthread \
          -Xcompiler -static-libstdc++ -Xcompiler -static-libgcc
else
GLUE_FLAGS =
GS_LIBS =
endif

.PHONY: all tools clean
all: gpusieve tools
tools: absort absplit ab2keys gsc2ab

# the glue's object differs with and without BOINC
GLUE_OBJ = boinc_glue$(if $(filter 1,$(HAVE_BOINC)),_boinc).o
$(GLUE_OBJ): boinc_glue.cpp boinc_glue.h Makefile
	$(CXX) -std=c++17 $(OPT) $(GLUE_FLAGS) -c -o $@ $<

gpusieve: gpusieve.cu gsc.h $(GLUE_OBJ) boinc_glue.h
	$(NVCC) $(NVCCFLAGS) -o $@ gpusieve.cu $(GLUE_OBJ) $(GS_LIBS)

absort: absort.cpp
	$(CXX) -std=c++17 $(OPT) -fopenmp -o $@ $<

absplit: absplit.c
	$(CC) $(OPT) -o $@ $<

gsc2ab: gsc2ab.cpp gsc.h
	$(CXX) -std=c++17 $(OPT) -o $@ $<

ab2keys: ab2keys.c
	$(CC) $(OPT) -o $@ $<

clean:
	rm -f gpusieve absort absplit ab2keys gsc2ab boinc_glue.o boinc_glue_boinc.o
