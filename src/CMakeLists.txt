aux_source_directory(. SRC_LIST)
# aux_source_directory(./Rmath RMATH)

include_directories(../include)
include_directories(../include/Rmath)

link_directories(../lib)

add_executable(main ${SRC_LIST})

target_link_libraries(main m Rmath)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
