# -*- mode: cmake -*-
#
# ParMETIS Find Module for Femus
#
# Search priority:
#
#   1. PARMETIS_INCLUDE_DIR / PARMETIS_LIBRARY_DIR
#   2. PETSc installation:
#        ${PETSC_DIR}/${PETSC_ARCH}/include
#        ${PETSC_DIR}/${PETSC_ARCH}/lib
#   3. PETSc externalpackages build/source directories
#   4. PARMETIS_DIR
#   5. Default CMake search paths
#
# Variables defined:
#
#   PARMETIS_FOUND
#   PARMETIS_INCLUDE_DIR
#   PARMETIS_INCLUDE_DIRS
#   PARMETIS_LIBRARY_DIR
#   PARMETIS_LIBRARY
#   PARMETIS_LIBRARIES
#

include(FindPackageHandleStandardArgs)


if(PARMETIS_LIBRARIES AND PARMETIS_INCLUDE_DIRS)

    # Already found. Do nothing.

else()

    # ------------------------------------------------------------------
    # Cache user-provided variables
    # ------------------------------------------------------------------

    if(PARMETIS_DIR)
        set(PARMETIS_DIR "${PARMETIS_DIR}"
            CACHE PATH "ParMETIS installation prefix")
    endif()

    if(PARMETIS_INCLUDE_DIR)
        set(PARMETIS_INCLUDE_DIR "${PARMETIS_INCLUDE_DIR}"
            CACHE PATH "ParMETIS include directory")
    endif()

    if(PARMETIS_LIBRARY_DIR)
        set(PARMETIS_LIBRARY_DIR "${PARMETIS_LIBRARY_DIR}"
            CACHE PATH "ParMETIS library directory")
    endif()


    # ==================================================================
    # INCLUDE
    # ==================================================================

    set(parmetis_inc_names "parmetis.h")


    # ------------------------------------------------------------------
    # First look where PETSc installs downloaded external packages.
    # This should normally be the correct location.
    # ------------------------------------------------------------------

    if(NOT PARMETIS_INCLUDE_DIR)

        find_path(
            PARMETIS_INCLUDE_DIR
            NAMES ${parmetis_inc_names}
            HINTS
                "${PETSC_DIR}/${PETSC_ARCH}/include"

                # Older / development PETSc externalpackage layouts
                "${PETSC_DIR}/${PETSC_ARCH}/externalpackages/git.parmetis/include"
                "${PETSC_DIR}/${PETSC_ARCH}/externalpackages/git.parmetis"
        )

    endif()


    # ------------------------------------------------------------------
    # If user explicitly provided PARMETIS_DIR
    # ------------------------------------------------------------------

    if(NOT PARMETIS_INCLUDE_DIR AND PARMETIS_DIR)

        find_path(
            PARMETIS_INCLUDE_DIR
            NAMES ${parmetis_inc_names}
            HINTS
                "${PARMETIS_DIR}"
            PATH_SUFFIXES
                include
                Include
            NO_DEFAULT_PATH
        )

    endif()


    # ------------------------------------------------------------------
    # Last fallback: normal CMake paths
    # ------------------------------------------------------------------

    if(NOT PARMETIS_INCLUDE_DIR)

        find_path(
            PARMETIS_INCLUDE_DIR
            NAMES ${parmetis_inc_names}
            PATH_SUFFIXES
                include
                Include
        )

    endif()


    if(NOT PARMETIS_INCLUDE_DIR)

        message(SEND_ERROR
            "Can not locate ParMETIS header parmetis.h. "
            "Checked PETSc at "
            "${PETSC_DIR}/${PETSC_ARCH}/include"
        )

    endif()


    # ==================================================================
    # LIBRARY
    # ==================================================================

    set(parmetis_lib_names "parmetis")


    # ------------------------------------------------------------------
    # PETSc normally installs libparmetis here.
    # ------------------------------------------------------------------

    if(NOT PARMETIS_LIBRARY_DIR)

        if(EXISTS "${PETSC_DIR}/${PETSC_ARCH}/lib")
            set(
                PARMETIS_LIBRARY_DIR
                "${PETSC_DIR}/${PETSC_ARCH}/lib"
            )
        endif()

    endif()


    # ------------------------------------------------------------------
    # Search explicitly in PARMETIS_LIBRARY_DIR first.
    # ------------------------------------------------------------------

    if(PARMETIS_LIBRARY_DIR)

        find_library(
            PARMETIS_LIBRARY
            NAMES ${parmetis_lib_names}
            HINTS
                "${PARMETIS_LIBRARY_DIR}"
            NO_DEFAULT_PATH
        )

    endif()


    # ------------------------------------------------------------------
    # Search PARMETIS_DIR if needed.
    # ------------------------------------------------------------------

    if(NOT PARMETIS_LIBRARY AND PARMETIS_DIR)

        find_library(
            PARMETIS_LIBRARY
            NAMES ${parmetis_lib_names}
            HINTS
                "${PARMETIS_DIR}"
            PATH_SUFFIXES
                lib
                lib64
                Lib
            NO_DEFAULT_PATH
        )

    endif()


    # ------------------------------------------------------------------
    # PETSc externalpackages fallback.
    #
    # Usually unnecessary because PETSc copies/installs the library into
    # ${PETSC_DIR}/${PETSC_ARCH}/lib.
    # ------------------------------------------------------------------

    if(NOT PARMETIS_LIBRARY)

        find_library(
            PARMETIS_LIBRARY
            NAMES ${parmetis_lib_names}
            HINTS
                "${PETSC_DIR}/${PETSC_ARCH}/externalpackages/git.parmetis/lib"
                "${PETSC_DIR}/${PETSC_ARCH}/externalpackages/git.parmetis/build/lib"
        )

    endif()


    # ------------------------------------------------------------------
    # Last fallback: standard CMake search
    # ------------------------------------------------------------------

    if(NOT PARMETIS_LIBRARY)

        find_library(
            PARMETIS_LIBRARY
            NAMES ${parmetis_lib_names}
        )

    endif()


    if(NOT PARMETIS_LIBRARY)

        message(SEND_ERROR
            "Can not locate ParMETIS library libparmetis. "
            "Checked PETSc at "
            "${PETSC_DIR}/${PETSC_ARCH}/lib"
        )

    endif()


    # ==================================================================
    # RESULT VARIABLES
    # ==================================================================

    set(
        PARMETIS_INCLUDE_DIRS
        ${PARMETIS_INCLUDE_DIR}
    )

    # ParMETIS requires METIS.
    #
    # Your project already has FindMETIS.cmake, therefore append
    # METIS_LIBRARIES here if it has already been found.
    if(METIS_LIBRARIES)

        set(
            PARMETIS_LIBRARIES
            ${PARMETIS_LIBRARY}
            ${METIS_LIBRARIES}
        )

    else()

        set(
            PARMETIS_LIBRARIES
            ${PARMETIS_LIBRARY}
        )

    endif()

endif()


# ======================================================================
# Standard CMake result
# ======================================================================

find_package_handle_standard_args(
    ParMETIS
    DEFAULT_MSG
    PARMETIS_LIBRARY
    PARMETIS_INCLUDE_DIR
)


# Femus convention uses uppercase PARMETIS_FOUND.
if(ParMETIS_FOUND)
    set(PARMETIS_FOUND TRUE)
else()
    set(PARMETIS_FOUND FALSE)
endif()


mark_as_advanced(
    PARMETIS_INCLUDE_DIR
    PARMETIS_INCLUDE_DIRS
    PARMETIS_LIBRARY
    PARMETIS_LIBRARIES
    PARMETIS_LIBRARY_DIR
)
