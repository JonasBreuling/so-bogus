
# required cmake version 
cmake_minimum_required(VERSION 2.8)

# packages
find_package(CUDA REQUIRED)

# nvcc flags
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_20,code=sm_20;-O6 )

include_directories( ../src )

cuda_compile(GS_O mat.cu)

cuda_add_executable(gs ${HELPER_O} ${GS_O} gs.cpp )
target_link_libraries( gs bogus )
target_link_libraries( gs ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcusparse.so )


