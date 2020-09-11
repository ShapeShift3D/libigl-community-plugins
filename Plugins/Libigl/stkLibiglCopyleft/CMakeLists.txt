if (NOT Eigen3_FOUND)
  vtk_module_find_package(PACKAGE Eigen3 VERSION 3.3 CONFIG_MODE PRIVATE)
endif ()

#MESSAGE ("LIBIGL_FOUND: " ${LIBIGL_FOUND})
if (NOT TARGET igl_common)
  vtk_module_find_package(PACKAGE LIBIGL COMPONENTS CORE TRIANGLE PRIVATE)
endif ()

if (NOT TARGET CGAL)
  vtk_module_find_package(PACKAGE CGAL PRIVATE)
endif ()

if(CGAL_FOUND)
  include( ${CGAL_USE_FILE} )
  
  # Activate concurrency ? (turned ON by default)
  option(CGAL_ACTIVATE_CONCURRENCY
         "Enable concurrency"
         ON)
   
  if(CGAL_ACTIVATE_CONCURRENCY)
    find_package(TBB REQUIRED)
    if(TBB_FOUND)
      add_definitions( -DNOMINMAX -DCGAL_LINKED_WITH_TBB )
    endif(TBB_FOUND)
  endif(CGAL_ACTIVATE_CONCURRENCY)
endif(CGAL_FOUND)

# ParaView 6.0 or later
set(classes
  stkLibiglBoolean3DMesher
  stkLibiglCopyLeftScalarClip
)

vtk_module_add_module(Libigl::stkLibiglCopyleft
  CLASSES ${classes})

paraview_add_server_manager_xmls(
  XMLS  stkLibiglCopyleftPlugins.xml)

vtk_module_include(Libigl::stkLibiglCopyleft PRIVATE
                    ${LIBIGL_INCLUDE_DIR})

vtk_module_link(Libigl::stkLibiglCopyleft PRIVATE CGAL::CGAL)

vtk_module_link(Libigl::stkLibiglCopyleft PUBLIC Eigen3::Eigen)

if(TBB_FOUND)
  message("TBB_FOUND: " ${TBB_FOUND})
  vtk_module_link(CGAL::vtkCGAL PRIVATE TBB::tbbmalloc TBB::tbb)
endif(TBB_FOUND)