cmake_minimum_required (VERSION 3.10)
project (mcts LANGUAGES CXX)

set(MCTS_BRIEF "Monte Carlo Tree Search framework")

if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "This application requires an out of source build.
        Please create a separate build directory")
endif()

add_library(${PROJECT_NAME} STATIC ${PROJECT_NAME}.hpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${PROJECT_NAME}.hpp")
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "${PROJECT_NAME}")
set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX ".a")
set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")

install( TARGETS "${PROJECT_NAME}"
   ARCHIVE DESTINATION "${INSTALL_DIR}/lib"
   COMPONENT library
   PUBLIC_HEADER DESTINATION "${INSTALL_DIR}/include"
)

target_compile_options    (${PROJECT_NAME} PRIVATE -O3 -Werror -Wall -Wextra -pedantic)
target_compile_features   (${PROJECT_NAME} PRIVATE cxx_std_17)