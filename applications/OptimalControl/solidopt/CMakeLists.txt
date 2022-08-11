CMAKE_MINIMUM_REQUIRED(VERSION 2.9)


# Full Dirichlet solid mechanics problem
ADD_SUBDIRECTORY(./solid_state) 

# distributed control solid mechanics problem
ADD_SUBDIRECTORY(./solid_control_lifting) 

# boundary control solid mechanics problem
ADD_SUBDIRECTORY(./solid_control_boundary) 
