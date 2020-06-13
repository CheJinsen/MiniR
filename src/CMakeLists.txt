aux_source_directory(. SRC_LIST)

include_directories(../include /usr/local/include)

link_directories(/usr/local/lib)

add_executable(stat ${SRC_LIST})

target_link_libraries(stat gsl gslcblas m)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
