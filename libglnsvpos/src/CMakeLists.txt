if (NOT DEFINED BUILD_SHARED_LIBS)
  set(BUILD_SHARED_LIBS ON)
endif (NOT DEFINED BUILD_SHARED_LIBS)

file(GLOB libglnsvpos_HEADERS "${PROJECT_SOURCE_DIR}/include/libglnsvpos/*.h")

include_directories("${PROJECT_SOURCE_DIR}/include")

set(libglnsvpos_SRCS
    glnsvpos.cpp
    rungekutta.cpp
)

add_library(glnsvpos-static STATIC ${libglnsvpos_SRCS})
install(TARGETS glnsvpos-static DESTINATION lib${LIB_SUFFIX})

if(BUILD_SHARED_LIBS)
  add_library(glnsvpos SHARED ${libglnsvpos_SRCS})
  install(TARGETS glnsvpos DESTINATION lib${LIB_SUFFIX})
else(BUILD_SHARED_LIBS)
  message(STATUS "Not building shared libraries")
endif(BUILD_SHARED_LIBS)

install(FILES ${libglnsvpos_HEADERS} DESTINATION include)
