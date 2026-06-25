
     set(VTUNE_HOME "/opt/intel/vtune_profiler")


     set(arch "64")

     find_path(VTUNE_INCLUDE ittnotify.h PATHS ${VTUNE_HOME}/include)
     find_library(VTUNE_LIBRARY libittnotify.a PATHS ${VTUNE_HOME}/lib${arch}/)

     if (NOT VTUNE_INCLUDE MATCHES NOTFOUND)
         if (NOT VTUNE_LIBRARY MATCHES NOTFOUND)
             set(VTUNE_FOUND TRUE)
             message(STATUS "ITT was found here ${VTUNE_HOME}")

             get_filename_component(VTUNE_LIBRARY_PATH ${VTUNE_LIBRARY} PATH)

             include_directories(${VTUNE_INCLUDE})
             link_directories(${VTUNE_LIBRARY_PATH})

             set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMFX_TRACE_ENABLE_ITT")
             set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMFX_TRACE_ENABLE_ITT")

             set(ITT_CFLAGS "-I${VTUNE_INCLUDE} -DITT_SUPPORT")
             set(ITT_LIBRARY_DIRS "${VTUNE_LIBRARY_PATH}")

             set(ITT_LIBS "")
             list(APPEND ITT_LIBS
                     mfx_trace
                     ittnotify+
                     )
         endif ()
     endif ()

