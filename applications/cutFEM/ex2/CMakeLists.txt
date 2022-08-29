CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

get_filename_component(APP_FOLDER_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

set( CMAKE_CXX_STANDARD 11 )

set(THIS_APPLICATION ${APP_FOLDER_NAME_PARENT}_${APP_FOLDER_NAME})

PROJECT(${THIS_APPLICATION})

ADD_EXECUTABLE(${THIS_APPLICATION}
${CMAKE_PROJECT_DIR}ex2.cpp 
${CMAKE_PROJECT_DIR}./../eqPoly/eqPoly.cpp
${CMAKE_PROJECT_DIR}./../eqPoly/bestFit.cpp
)


find_library( CLN_LIB "cln" DOC "Path to libcln" )
target_link_libraries(${THIS_APPLICATION} ${CLN_LIB} )

TARGET_LINK_LIBRARIES(${THIS_APPLICATION} femus)
TARGET_LINK_LIBRARIES(${THIS_APPLICATION} ${PETSC_LIBRARIES})
TARGET_LINK_LIBRARIES(${THIS_APPLICATION} b64)
TARGET_LINK_LIBRARIES(${THIS_APPLICATION} jsoncpp)
TARGET_LINK_LIBRARIES(${THIS_APPLICATION} adept)

IF(FPARSER_FOUND)
  TARGET_LINK_LIBRARIES(${THIS_APPLICATION} ${FPARSER_LIBRARY})
ENDIF(FPARSER_FOUND)

IF(MPI_FOUND) 
  TARGET_LINK_LIBRARIES(${THIS_APPLICATION} ${MPI_CXX_LIBRARIES})
ENDIF(MPI_FOUND)

IF(SLEPC_FOUND)
  TARGET_LINK_LIBRARIES(${THIS_APPLICATION}  ${SLEPC_LIBARIES})
  TARGET_LINK_LIBRARIES(${THIS_APPLICATION}  SLEPC::slepc SLEPC::slepc_static)
ENDIF(SLEPC_FOUND)

IF(HDF5_FOUND)
  TARGET_LINK_LIBRARIES(${appname} ${HDF5_LIBRARIES})
ENDIF(HDF5_FOUND)

FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output/)
FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/input/)
FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/save/)
FILE(COPY           ${PROJECT_SOURCE_DIR}/input/ DESTINATION ${PROJECT_BINARY_DIR}/input/)

# ### Minimum required version
# cmake_minimum_required( VERSION 3.3 )
# 
# ### Project name
# project( "LiSK" VERSION 1.0 LANGUAGES CXX )
# 
# ### Special commands for compilation
# set( CMAKE_CXX_STANDARD 11 )
# 
# ### Add CLN library
# ### Header file and library should automatically be
# ### found if present in environment
# find_path( CLN_INCLUDE_DIR "cln/cln.h" DOC "Directory of cln.h" )
# find_library( CLN_LIB "cln" DOC "Path to libcln" )
# 
# ### Executable
# add_executable( example "eqPoly.cpp" )
# set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
# 
# ### Include directories
# include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../LiSK )
# include_directories( SYSTEM ${CLN_INCLUDE_DIR} )
# 
# ### Linking
# target_link_libraries( example ${CLN_LIB} )
