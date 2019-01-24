# ======================================================================
#  larcvconversion main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larcvconversion/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================


# use cmake 2.8 or later
cmake_minimum_required (VERSION 2.8)

project(larcvconversion)

# cetbuildtools contains our cmake modules
SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
IF (NOT CETBUILDTOOLS_VERSION)
    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules")
ENDIF()

set(CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules
		      $ENV{CETBUILDTOOLS_DIR}/Modules
		      ${CMAKE_MODULE_PATH})

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  ALLOW_DEPRECATIONS
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs  -O0
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_product( larsoft v1_00_00 )
find_ups_product( nutools v1_00_00 )
find_ups_product( art v1_08_00 )
find_ups_product( cetbuildtools v4_14_01 )

#  Find all the libraries needed by our dependent CMakeList.txt files
cet_find_library( BOOST_SERIALIZATION NAMES boost_serialization PATHS ENV BOOST_LIB NO_DEFAULT_PATH )
cet_find_library( BOOST_DATE_TIME     NAMES boost_date_time     PATHS ENV BOOST_LIB NO_DEFAULT_PATH )


if(DEFINED ENV{LAPACK_LIB} )
  find_ups_product( lapack )
endif()

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# ADD SOURCE CODE SUBDIRECTORIES HERE
add_subdirectory(larcvconversion)

### tests
###add_subdirectory(test)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)
