if (NOT Eigen3_FOUND)
  vtk_module_find_package(PACKAGE Eigen3 VERSION 3.3 CONFIG_MODE PRIVATE)
endif ()

if (${LIBIGL_DIR} STREQUAL "")
  message(FATAL_ERROR "Path to LIBIGL_DIR is empty")
endif()

vtk_module_find_package(PACKAGE CGAL VERSION 5.1.1 CONFIG_MODE PRIVATE)

# Namespace adaptive to how the module is built (as a pure module or as plugins)
SET(NAMESPACE "")

if(ENABLE_NONFREE_MODULES)
  SET(NAMESPACE "ParaView")
endif(ENABLE_NONFREE_MODULES)
  
if(ENABLE_NONFREE_PLUGINS)
  SET(NAMESPACE "Libigl")
endif(ENABLE_NONFREE_PLUGINS)

# Define content
set(classes
  stkLibiglBoolean3DMesher
)
  
vtk_module_add_module(${NAMESPACE}::stkLibiglCopyleft
  CLASSES ${classes}
  PRIVATE_HEADERS ${private_headers}
)

set(xml_files
  Resources/LibiglCopyleftProxies.xml
  Resources/LibiglCopyleftFilters.xml
)

paraview_add_server_manager_xmls(
  XMLS  ${xml_files}
)

if(CGAL_FOUND)
  include( ${CGAL_USE_FILE} )
  
  # Activate concurrency ? (turned ON by default)
  option(CGAL_ACTIVATE_CONCURRENCY
         "Enable concurrency"
         ON)
   
  vtk_module_definitions(${NAMESPACE}::stkLibiglCopyleft PRIVATE NOMINMAX)

  # Link with CGAL
  vtk_module_link(${NAMESPACE}::stkLibiglCopyleft PUBLIC CGAL::CGAL)

  if(CGAL_ACTIVATE_CONCURRENCY)
    find_package(TBB REQUIRED)
    if(TBB_FOUND)
      vtk_module_definitions(${NAMESPACE}::stkLibiglCopyleft PRIVATE CGAL_LINKED_WITH_TBB)
      vtk_module_link(${NAMESPACE}::stkLibiglCopyleft PRIVATE TBB::tbbmalloc TBB::tbb)
    endif(TBB_FOUND)
  endif(CGAL_ACTIVATE_CONCURRENCY)
endif(CGAL_FOUND)

vtk_module_include(${NAMESPACE}::stkLibiglCopyleft PRIVATE ${LIBIGL_DIR}/include)

vtk_module_link(${NAMESPACE}::stkLibiglCopyleft PUBLIC Eigen3::Eigen)