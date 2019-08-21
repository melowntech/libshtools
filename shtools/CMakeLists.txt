# bump version here
set(shtools_VERSION 1.0)

define_module(LIBRARY shtools=${shtools_VERSION}
  DEPENDS BuildSystem>=1.6 Fortran C++>=14 OpenMP Boost)

# grap all sources sans python extension
file(GLOB shtools_f95_SOURCES CONFIGURE_DEPENDS shtools/*.f95 shtools/*.F95)

list(LENGTH shtools_f95_SOURCES shtools_f95_SOURCES_LENGTH)
if(shtools_f95_SOURCES_LENGTH EQUAL 0)
  message(fatal "Missing SHTOOLS sources. Did you foreget to check"
    "out SHTOOLS submodule?")
endif()

set(shtools_SOURCES
  shtools.hpp shtools.cpp
  )

add_library(shtools STATIC ${shtools_SOURCES} ${shtools_f95_SOURCES})
buildsys_library(shtools)
# TODO: use cmake machinery around those libraries
target_link_libraries(shtools ${MODULE_LIBRARIES} -lfftw3 -llapack -lblas)
buildsys_target_compile_definitions(shtools ${MODULE_DEFINITIONS})

