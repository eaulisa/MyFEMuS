CMAKE_MINIMUM_REQUIRED(VERSION 3.15)

get_filename_component(APP_FOLDER_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

set(THIS_APPLICATION ${APP_FOLDER_NAME})

PROJECT("${APP_FOLDER_NAME_PARENT}_${THIS_APPLICATION}")


SET(MAIN_FILE "${THIS_APPLICATION}") # the name of the main file with no extension
SET(EXEC_FILE "${APP_FOLDER_NAME_PARENT}_${MAIN_FILE}") # the name of the executable file

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

femusMacroBuildApplication(${MAIN_FILE} ${EXEC_FILE})

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_path(GPERFTOOLS_INCLUDE_DIR profiler.h PATH_SUFFIXES gperftools)
find_library(GPERFTOOLS_LIB profiler)

if (GPERFTOOLS_INCLUDE_DIR AND GPERFTOOLS_LIB)
    message(STATUS "Found gperftools: ${GPERFTOOLS_LIB}")
    include_directories(${GPERFTOOLS_INCLUDE_DIR})
    target_link_libraries(${EXEC_FILE} ${GPERFTOOLS_LIB})
else()
    message(WARNING "gperftools not found, profiling disabled")
endif()

FILE(COPY ${PROJECT_SOURCE_DIR}/paraviewReadMarker.xml DESTINATION ${PROJECT_BINARY_DIR}/)
