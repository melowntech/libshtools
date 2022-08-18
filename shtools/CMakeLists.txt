# bump version here
set(shtools_VERSION 1.1)

if(BUILDSYS_CONAN)
  list(APPEND shtools_EXTRA_MODULE_LIBRARIES shtools::shtools)
else()
  list(APPEND shtools_EXTRA_DEPENDS Fortran)
endif()

define_module(LIBRARY shtools=${shtools_VERSION}
  DEPENDS BuildSystem>=1.13 C++>=14 OpenMP Boost
  LAPACK BLAS FFTW3 ${shtools_EXTRA_DEPENDS})

if(NOT BUILDSYS_CONAN)
  # grap all sources sans python extension
  file(GLOB shtools_f95_SOURCES CONFIGURE_DEPENDS shtools/*.f95 shtools/*.F95)

  list(LENGTH shtools_f95_SOURCES shtools_f95_SOURCES_LENGTH)
  if(shtools_f95_SOURCES_LENGTH EQUAL 0)
    message(FATAL_ERROR "Missing SHTOOLS sources. Did you foreget to check"
      "out SHTOOLS submodule?")
  endif()
endif()

set(shtools_SOURCES
  shtools.hpp shtools.cpp
  )

add_library(shtools STATIC ${shtools_SOURCES} ${shtools_f95_SOURCES})
buildsys_library(shtools)
# TODO: use cmake machinery around those libraries
target_link_libraries(shtools ${MODULE_LIBRARIES} ${shtools_EXTRA_MODULE_LIBRARIES})
target_compile_definitions(shtools PRIVATE ${MODULE_DEFINITIONS})

