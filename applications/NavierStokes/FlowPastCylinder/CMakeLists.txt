CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(TimedepNavierStokes)

SET(MAIN_FILE "main")
SET(EXEC_FILE "flowPastCylinder")

femusMacroBuildApplication(${MAIN_FILE} ${EXEC_FILE})

