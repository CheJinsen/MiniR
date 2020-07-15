aux_source_directory(. SRC_LIST)

include_directories(../include)

# link_directories(/usr/local/lib)

add_executable(main ${SRC_LIST})

target_link_libraries(main m stdc++)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
