CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(TimedepNavierStokes)

SET(MAIN_FILE "main")
SET(EXEC_FILE "flowPastCylinder")

femusMacroBuildApplication(${MAIN_FILE} ${EXEC_FILE})

