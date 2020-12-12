aux_source_directory(. SRC_LIST)
aux_source_directory(./Rmath RMATH)

include_directories(../include)
include_directories(../include/Rmath)

# link_directories(/usr/local/lib)

add_executable(main ${SRC_LIST} ${RMATH})

target_link_libraries(main m)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
