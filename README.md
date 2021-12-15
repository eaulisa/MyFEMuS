MyFEMuS
======

Welcome to the MyFEMuS project! MyFEMuS is a fork of the FEMuS project administered mainly by Eugenio Aulisa.
The manual intallation, will also work with the FEMuS project.

For the FEMuS project automatic installation, as well as the Mac-specific installation, see below.


<!-- ![alt tag](https://github.com/FeMTTU/femus/blob/master/doc/images/logo.jpg?raw=true) -->
<!-- ![alt tag](https://github.com/FeMTTU/femus/blob/master/doc/images/FSI.jpg?raw=true) -->

Step by step MyFEMuS manual setup  (largely tested on OpenSuse)
======

Install PETSC

From the directory $INSTALLATION_DIR clone petsc

    git clone -b maint https://bitbucket.org/petsc/petsc petsc
    
    cd petsc
    
Configure, compile, test PETSC with the following options
    
    ./configure --with-debugging=0 --with-x=1 COPTFLAGS="-O3 -march=native -mtune=native" CXXOPTFLAGS="-O3 -march=native -mtune=native" FOPTFLAGS="-O3 -march=native -mtune=native" --download-openmpi=1 --download-fblaslapack=1 --download-hdf5=1 --download-metis=1 --download-parmetis=1 --with-shared-libraries=1 --download-blacs=1 --download-scalapack=1 --download-mumps=1 --download-suitesparse

    make PETSC_DIR=$INSTALLATION_DIR/petsc PETSC_ARCH=arch-linux2-c-opt all

    make PETSC_DIR=$INSTALLATION_DIR/petsc PETSC_ARCH=arch-linux2-c-opt test
 
======

Install SLEPC

From the directory $INSTALLATION_DIR clone slepc

    git clone -b maint https://bitbucket.org/slepc/slepc slepc

    cd slepc
    
Configure, compile, test SLEPC with the following options
    
    export PETSC_DIR=$INSTALLATION_DIR/petsc 
    
    export PETSC_ARCH=arch-linux2-c-opt
    
    ./configure
    
    make SLEPC_DIR=$PWD all
    
    make SLEPC_DIR=$PWD test

======

Install MyFEMuS 

Be sure you have installed al least gcc 7, cmake, cmake-gui. Fparser may be handy for some applications but it is not required. 

Clone the MyFEMuS source code from the github repository

From the directory $INSTALLATION_DIR clone MyFEMuS

    https://github.com/eaulisa/MyFEMuS.git
    
    cd MyFEMuS
    
I generally export the following variables in the ./bashrc file in my user home, so that are available everywhere, otherwise you will need to export them all the times. 
    
    export PETSC_DIR=$INSTALLATION_DIR/petsc 
    
    export PETSC_ARCH=arch-linux2-c-opt
    
    export SLEPC_DIR=$INSTALLATION_DIR/slepc 
    
Configure MyFEMuS using cmake-gui. 

    cmake-gui 

    Where is the source code: $INSTALLATION_DIR/MyFEMuS
    
    Where to build the binaries: $INSTALLATION_DIR/feumsbin
    
    CMAKE_BUILD_TYPE choose between release (default) or debug
    
    Press Configure button
    
    Press Generate button

Compile
    
    cd $INSTALLATION_DIR/femusbin
    
    make
    
Run. All applications are built in the folder $INSTALLATION_DIR/femusbin/applications/..
    
======
    

FEMuS automatic configuration, contact Giorgio Bornia for support.
======

Welcome to the FEMuS project! FEMuS is an open-source Finite Element C++ library 
built on top of PETSc, which allows scientists to build and solve multiphysics 
problems with multigrid and domain decomposition techniques.


<!-- ![alt tag](https://github.com/FeMTTU/femus/blob/master/doc/images/logo.jpg?raw=true) -->
<!-- ![alt tag](https://github.com/FeMTTU/femus/blob/master/doc/images/FSI.jpg?raw=true) -->

Setup
=====


Clone the FEMuS source code from the github repository:


    git clone https://github.com/FeMTTU/femus.git

   
You need PETSc for FEMuS to work.
If PETSc is not already installed in your machine, the script "install_petsc.sh" in contrib/scripts/ will install it automatically,
with the following syntax:
  
    ./femus/contrib/scripts/install_petsc.sh --prefix-external my_dir 
  

where "my_dir" is the directory, either absolute or relative, in which you want PETSc to be installed 
(please put it outside of the femus repo directory, to prevent from potential git tracking).

 Source the "configure_femus.sh" script and execute the function "fm_set_femus" in order to set some environment variables:

    source femus/contrib/scripts/configure_femus.sh

    fm_set_femus  --prefix-external my_dir --method-petsc opt
   

  Create the build directory, cd to it and run cmake:
   
    mkdir femus.build

    cd femus.build

    cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE="[Debug Release RelWithDebInfo MinSizeRel None]"  ../femus

=====


FEMuS Mac installation, contact Anthony Gruber for support.  Note that the optional FParser and Libmesh functionality is not yet usable.
======

Download homebrew

    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

Install the following packages 

    Brew install gcc
    Brew install make
    Brew install open-mpi
    Brew install metis
    Brew install parmetis
    Brew install hdf5
    Brew install scalapack
    Brew install boost
    Brew install --cask cmake (for cmake-gui)
    Brew install pkg-config (maybe not necessary)

Install PETSc

From the directory $INSTALLATION_DIR clone petsc

    git clone -b release https://gitlab.com/petsc/petsc.git petsc 
    
    cd petsc
    
Configure PETSc with the following options (tested 12/14/2021 on MBP 2021 -- replace directories according to your own homebrew installations)
    
    ./configure --with-debugging=0 --with-shared-libraries --with-mpi-dir=/opt/homebrew/Cellar/open-mpi/4.1.2 --with-hdf5-dir=/opt/homebrew/Cellar/hdf5/1.12.1 --with-boost-dir=/opt/homebrew/Cellar/boost/1.76.0 --with-metis-dir=/opt/homebrew/Cellar/metis/5.1.0 --with-parmetis-dir=/opt/homebrew/Cellar/parmetis/4.0.3_5 --with-scalapack-dir=/opt/homebrew/Cellar/scalapack/2.1.0_3 --download-mumps --download-blacs --download-suitesparse

Follow the console prompts to compile and test the PETSc library.  All tests should pass.

Install SLEPc

From the directory $INSTALLATION_DIR clone slepc

    git clone -b release https://gitlab.com/slepc/slepc

Put the following lines in your .zshrc (or just run them in local scope from the shell)

    export PETSC_DIR=$INSTALLATION_DIR/petsc 

    export PETSC_ARCH=arch-darwin-c-opt

    export SLEPC_DIR=$INSTALLATION_DIR/slepc

Configure SLEPc

    ./configure

Follow the console prompts to compile and test the SLEPc library.  All tests should pass.

Install MyFEMuS

From the directory $INSTALLATION_DIR clone MyFEMuS and make a directory for the binaries

    git clone https://github.com/agrubertx/MyFEMuS.git

    mkdir femusbin

Navigate to MyFEMuS and checkout branch "anthony"

    cd MyFEMuS

    git checkout anthony

Configure MyFEMuS using cmake-gui. 

    cmake-gui 

    Where is the source code: $INSTALLATION_DIR/MyFEMuS
    
    Where to build the binaries: $INSTALLATION_DIR/feumsbin
    
    CMAKE_BUILD_TYPE choose between release (default) or debug
    
    Press Configure button
    
    Press Generate button

Compile
    
    cd $INSTALLATION_DIR/femusbin
    
    make
    
Run. All applications are built in the folder $INSTALLATION_DIR/femusbin/applications/..



Authors
========

Eugenio Aulisa

Simone Bnà

Giorgio Bornia

Anthony Gruber



License
========

FEMuS is an open-source software distributed under the LGPL license, version 2.1

