# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/pacMapper/CMakeLists.txt
#
# CMakeLists.txt file for INDEX_BENCH.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_apps_pacMapper)
message (STATUS "Configuring apps/qbin")

set (SEQAN_APP_VERSION "1.1.1")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------
find_package(CXX11 REQUIRED)
if (NOT CXX11_FOUND)
    message (STATUS " C++11 compiler needed to build qbin")
    return ()
endif (NOT CXX11_FOUND)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message (STATUS "OpenMP needed to build qbin")
endif()

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# Update the list of file names below if you add source files to your application.
add_executable (qbin testBinning.cpp mapper.h base.h index_extend.h shape_extend.h mapparm.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (qbin ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${CXX11_CXX_FLAGS}")
message ("debug cxx" ${CMAKE_CXX_FLAGS})

set (CMAKE_BUILD_TYPE Release)

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
#if (NOT SEQAN_PREFIX_SHARE_DOC)
#  seqan_setup_install_vars (qbin)
#endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install pacMapper in ${PREFIX}/bin directory
install (TARGETS qbin
         DESTINATION bin/)
# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/pacMapper
#install (FILES LICENSE
#               README
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/small.fasta
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (pacMapper)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable seqan_tcoffee in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} qbin CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:qbin")
  set (CPACK_PACKAGE_NAME "qbin")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "qbin - qindex for binning")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Chenxu Pan <chenxu.pan@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (qbin "qbin")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:qbin")
