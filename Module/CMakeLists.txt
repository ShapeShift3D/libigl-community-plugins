if(NOT Eigen3_FOUND)
  vtk_module_find_package(PACKAGE Eigen3 VERSION 3.3 CONFIG_MODE PRIVATE)
endif()

vtk_module_find_package(PACKAGE CGAL VERSION 5.1.1 CONFIG_MODE PRIVATE)

if(NOT DEFINED LIBIGL_DIR)
  message(FATAL_ERROR "The CMake variable LIBIGL_DIR is undefined.")
endif()

if(${LIBIGL_DIR} STREQUAL "")
  message(FATAL_ERROR "The CMake variable LIBIGL_DIR is an empty string.")
endif()

if(NOT TARGET igl::core)
  include(LibiglOptions)
  add_subdirectory(${LIBIGL_DIR} ${CMAKE_BINARY_DIR}/bin/libigl)
endif()

# Define content
set(classes
  stkLibiglBoolean3DMesher
)

vtk_module_add_module(STK::LibiglCopyleft
  CLASSES ${classes}
  PRIVATE_HEADERS ${private_headers}
)

if(CGAL_FOUND)
  include(${CGAL_USE_FILE})

  # Activate concurrency ? (turned ON by default)
  option(CGAL_ACTIVATE_CONCURRENCY
    "Enable concurrency"
    ON)

  vtk_module_definitions(STK::LibiglCopyleft PRIVATE NOMINMAX)

  # Link with CGAL
  vtk_module_link(STK::LibiglCopyleft PUBLIC CGAL::CGAL)

  if(CGAL_ACTIVATE_CONCURRENCY)
    if (NOT TARGET TBB::tbb)
      find_package(TBB REQUIRED)
    endif ()

    if(TBB_FOUND)
      vtk_module_definitions(STK::LibiglCopyleft PRIVATE CGAL_LINKED_WITH_TBB)
      vtk_module_link(STK::LibiglCopyleft PRIVATE TBB::tbbmalloc TBB::tbb)
    endif(TBB_FOUND)
  endif(CGAL_ACTIVATE_CONCURRENCY)
endif(CGAL_FOUND)

vtk_module_include(STK::LibiglCopyleft PRIVATE ${LIBIGL_DIR}/include)

vtk_module_link(STK::LibiglCopyleft PUBLIC Eigen3::Eigen)