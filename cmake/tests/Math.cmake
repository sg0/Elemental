# Check for BLAS and LAPACK support
# =================================
if(EL_PURE)
  set(MATH_DESC "Unthreaded BLAS/LAPACK link flags")
else()
  set(MATH_DESC "Threaded BLAS/LAPACK link flags")
endif()
if(MATH_LIBS)
  message(STATUS "Using user-defined MATH_LIBS=${MATH_LIBS}")
elseif(APPLE)
  # Try to determine whether to use vecLib (old) or Accelerate (new)
  set(CMAKE_REQUIRED_FLAGS "${MPI_C_COMPILE_FLAGS} ${MPI_LINK_FLAGS}")
  set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
  set(CMAKE_REQUIRED_LIBRARIES "-framework vecLib;${MPI_C_LIBRARIES}")
  check_function_exists(dpotrf  EL_HAVE_DPOTRF_VECLIB)
  check_function_exists(dpotrf_ EL_HAVE_DPOTRF_POST_VECLIB)
  set(CMAKE_REQUIRED_LIBRARIES "-framework Accelerate;${MPI_C_LIBRARIES}")
  check_function_exists(dpotrf  EL_HAVE_DPOTRF_ACCELERATE)
  check_function_exists(dpotrf_ EL_HAVE_DPOTRF_POST_ACCELERATE)
  if(EL_HAVE_DPOTRF_VECLIB OR EL_HAVE_DPOTRF_POST_VECLIB)
    set(MATH_LIBS "-framework vecLib" CACHE STRING ${MATH_DESC})
    message(STATUS "Using Apple vecLib framework.")
  elseif(EL_HAVE_DPOTRF_ACCELERATE OR EL_HAVE_DPOTRF_POST_ACCELERATE)
    set(MATH_LIBS "-framework Accelerate" CACHE STRING ${MATH_DESC})
    message(STATUS "Using Apple Accelerate framework.")
  endif()
else()
  # Look for default BLAS and LAPACK
  if(REFERENCE_ROOT)
    message(STATUS "Searching REFERENCE_ROOT=${REFERENCE_ROOT} for math libs")
  endif()
  set(REFERENCE_REQUIRED LAPACK BLAS)
  find_library(BLAS_LIB NAMES blas PATHS ${REFERENCE_ROOT})
  find_library(LAPACK_LIB NAMES lapack reflapack PATHS ${REFERENCE_ROOT})
  set(REFERENCE_FOUND TRUE)
  foreach(NAME ${REFERENCE_REQUIRED})
    if(${NAME}_LIB)
      message(STATUS "Found ${NAME}_LIB: ${${NAME}_LIB}")
      list(APPEND MATH_LIBS ${${NAME}_LIB})
    else()
      message(STATUS "Could not find ${NAME}_LIB")
      set(MATH_LIBS "")
      set(REFERENCE_FOUND FALSE)
    endif()
  endforeach()
  if(REFERENCE_FOUND)
    message(WARNING "Using reference BLAS/LAPACK; performance will be poor")
  else()
    message(FATAL_ERROR "Could not find BLAS/LAPACK. Please specify MATH_LIBS")
  endif()
endif()

# Check the BLAS and LAPACK underscore conventions
set(CMAKE_REQUIRED_FLAGS "${MPI_C_COMPILE_FLAGS} ${MPI_LINK_FLAGS}")
set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
set(CMAKE_REQUIRED_LIBRARIES "${MATH_LIBS};${MPI_C_LIBRARIES}")
check_function_exists(daxpy  EL_HAVE_DAXPY)
check_function_exists(daxpy_ EL_HAVE_DAXPY_POST)
check_function_exists(dpotrf  EL_HAVE_DPOTRF)
check_function_exists(dpotrf_ EL_HAVE_DPOTRF_POST)
if(EL_HAVE_DAXPY)
  set(EL_BLAS_POST FALSE)
  set(EL_BLAS_DEFS "")
elseif(EL_HAVE_DAXPY_POST)
  set(EL_BLAS_POST TRUE)
  set(EL_BLAS_DEFS "-DBLAS_POST")
else()
  message(FATAL_ERROR "Could not determine BLAS format.")
endif()
if(EL_HAVE_DPOTRF)
  set(EL_LAPACK_POST FALSE)
  set(EL_LAPACK_DEFS "")
elseif(EL_HAVE_DPOTRF_POST)
  set(EL_LAPACK_POST TRUE)
  set(EL_LAPACK_DEFS "-DLAPACK_POST")
else()
  message(FATAL_ERROR "Could not determine LAPACK format.")
endif()
# Ensure that we have a relatively new version of LAPACK
if(EL_HAVE_DPOTRF)
  check_function_exists(dsyevr EL_HAVE_DSYEVR)
  if(NOT EL_HAVE_DSYEVR)
    message(FATAL_ERROR "LAPACK is missing dsyevr")
  endif()
else()
  check_function_exists(dsyevr_ EL_HAVE_DSYEVR_POST)
  if(NOT EL_HAVE_DSYEVR_POST)
    message(FATAL_ERROR "LAPACK is missing dsyevr_")
  endif()
endif()

# Check for libFLAME support
# ==========================
check_function_exists(FLA_Bsvd_v_opd_var1 EL_HAVE_FLA_BSVD)

# Check for ScaLAPACK support
# ===========================
if(NOT EL_DISABLE_SCALAPACK)
  # NOTE: pdsyngst was chosen because MKL's ScaLAPACK only defines pdsyngst_,
  #       but not pdsyngst, despite defining both pdpotrf and pdpotrf_. 
  check_function_exists(pdsyngst  EL_HAVE_PDSYNGST)
  check_function_exists(pdsyngst_ EL_HAVE_PDSYNGST_POST)
  check_function_exists(Csys2blacs_handle EL_HAVE_CSYS2BLACS)
  if(EL_HAVE_PDSYNGST)
    check_function_exists(pdlaqr0 EL_HAVE_PDLAQR0)
    check_function_exists(pdlaqr1 EL_HAVE_PDLAQR1)
    if(NOT EL_HAVE_PDLAQR0 OR NOT EL_HAVE_PDLAQR1 OR 
       NOT EL_HAVE_CSYS2BLACS)
      message(STATUS "ScaLAPACK must support PDLAQR{0,1} and Csys2blacs_handle.")
      set(EL_HAVE_SCALAPACK FALSE)
    else()
      set(EL_HAVE_SCALAPACK TRUE)
    endif()
    set(EL_SCALAPACK_POST FALSE)
    set(EL_SCALAPACK_DEFS "")
  elseif(EL_HAVE_PDSYNGST_POST)
    check_function_exists(pdlaqr0_ EL_HAVE_PDLAQR0_POST)
    check_function_exists(pdlaqr1_ EL_HAVE_PDLAQR1_POST)
    if(NOT EL_HAVE_PDLAQR0_POST OR NOT EL_HAVE_PDLAQR1_POST OR
       NOT EL_HAVE_CSYS2BLACS)
      message(STATUS "ScaLAPACK must support PDLAQR{0,1} and Csys2blacs_handle.")
    else()
      set(EL_HAVE_SCALAPACK TRUE)
    endif()
    set(EL_SCALAPACK_POST TRUE)
    set(EL_SCALAPACK_DEFS "-DBLAS_POST")
  else()
    set(EL_HAVE_SCALAPACK FALSE)
    message(STATUS "ScaLAPACK was NOT detected.")
  endif()
endif()
