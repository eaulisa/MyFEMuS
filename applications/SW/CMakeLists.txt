CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

get_filename_component(APP_FOLDER_NAME_PARENT ${CMAKE_CURRENT_LIST_DIR} NAME)

IF(SLEPC_FOUND)

  ADD_SUBDIRECTORY(ex1/)

  ADD_SUBDIRECTORY(ex2/)

  ADD_SUBDIRECTORY(ex3/)

  ADD_SUBDIRECTORY(ex4/)

  ADD_SUBDIRECTORY(ex5/)

  ADD_SUBDIRECTORY(lock_exchange_isopycnal/)

  ADD_SUBDIRECTORY(lock_exchange_zlevel/)

  ADD_SUBDIRECTORY(lock_exchange_zlevel_split/)

  ADD_SUBDIRECTORY(lock_exchange_zlevel_2matrices/)

  ADD_SUBDIRECTORY(overflow_isopycnal/)

  ADD_SUBDIRECTORY(internal_waves_zlevel/)

  ADD_SUBDIRECTORY(overflow_zlevel/)

  ADD_SUBDIRECTORY(tracer/)
ENDIF(SLEPC_FOUND)
