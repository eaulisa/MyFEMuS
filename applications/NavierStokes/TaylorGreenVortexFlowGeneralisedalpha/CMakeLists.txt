CMAKE_MINIMUM_REQUIRED(VERSION 3.15)

PROJECT(TaylorVortex)

SET(MAIN_FILE "main")
SET(EXEC_FILE "taylorvortex")

femusMacroBuildApplication(${MAIN_FILE} ${EXEC_FILE})

