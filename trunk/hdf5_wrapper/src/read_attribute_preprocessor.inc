!
! Routines to read n-dimensional integer*4 datasets

#define NATIVE_TYPE NATIVE_INTEGER4
#define ARR_TYPE INTEGER*4
#define DATA_CLASS H5T_INTEGER_F

#define SUB_NAME read_1d_integer4_array_attribute
#define NDIMS 1
#define ARRAY_DIM (:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_2d_integer4_array_attribute
#define NDIMS 2
#define ARRAY_DIM (:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_3d_integer4_array_attribute
#define NDIMS 3
#define ARRAY_DIM (:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_4d_integer4_array_attribute
#define NDIMS 4
#define ARRAY_DIM (:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_5d_integer4_array_attribute
#define NDIMS 5
#define ARRAY_DIM (:,:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_6d_integer4_array_attribute
#define NDIMS 6
#define ARRAY_DIM (:,:,:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_7d_integer4_array_attribute
#define NDIMS 7
#define ARRAY_DIM (:,:,:,:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_scalar_integer4_attribute
#define NDIMS 0
#define SCALAR 1
#include "read_attribute.inc"
#undef SCALAR
#undef NDIMS
#undef SUB_NAME

#undef NATIVE_TYPE
#undef ARR_TYPE
#undef DATA_CLASS

!
! Routines to read n-dimensional integer*8 datasets

#define NATIVE_TYPE NATIVE_INTEGER8
#define ARR_TYPE INTEGER*8
#define DATA_CLASS H5T_INTEGER_F

#define SUB_NAME read_1d_integer8_array_attribute
#define NDIMS 1
#define ARRAY_DIM (:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_2d_integer8_array_attribute
#define NDIMS 2
#define ARRAY_DIM (:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_3d_integer8_array_attribute
#define NDIMS 3
#define ARRAY_DIM (:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_4d_integer8_array_attribute
#define NDIMS 4
#define ARRAY_DIM (:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_5d_integer8_array_attribute
#define NDIMS 5
#define ARRAY_DIM (:,:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_6d_integer8_array_attribute
#define NDIMS 6
#define ARRAY_DIM (:,:,:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_7d_integer8_array_attribute
#define NDIMS 7
#define ARRAY_DIM (:,:,:,:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_scalar_integer8_attribute
#define NDIMS 0
#define SCALAR 1
#include "read_attribute.inc"
#undef SCALAR
#undef NDIMS
#undef SUB_NAME

#undef NATIVE_TYPE
#undef ARR_TYPE
#undef DATA_CLASS

!
! Routines to read n-dimensional real datasets
#define NATIVE_TYPE NATIVE_REAL4
#define ARR_TYPE REAL*4
#define DATA_CLASS H5T_FLOAT_F

#define SUB_NAME read_1d_real_array_attribute
#define NDIMS 1
#define ARRAY_DIM (:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_2d_real_array_attribute
#define NDIMS 2
#define ARRAY_DIM (:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_3d_real_array_attribute
#define NDIMS 3
#define ARRAY_DIM (:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_4d_real_array_attribute
#define NDIMS 4
#define ARRAY_DIM (:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_5d_real_array_attribute
#define NDIMS 5
#define ARRAY_DIM (:,:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_6d_real_array_attribute
#define NDIMS 6
#define ARRAY_DIM (:,:,:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_7d_real_array_attribute
#define NDIMS 7
#define ARRAY_DIM (:,:,:,:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_scalar_real_attribute
#define NDIMS 0
#define SCALAR 1
#include "read_attribute.inc"
#undef SCALAR
#undef NDIMS
#undef SUB_NAME

#undef NATIVE_TYPE
#undef ARR_TYPE
#undef DATA_CLASS

!
! Routines to read n-dimensional real datasets
#define NATIVE_TYPE NATIVE_REAL8
#define ARR_TYPE REAL*8
#define DATA_CLASS H5T_FLOAT_F

#define SUB_NAME read_1d_double_array_attribute
#define NDIMS 1
#define ARRAY_DIM (:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_2d_double_array_attribute
#define NDIMS 2
#define ARRAY_DIM (:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_3d_double_array_attribute
#define NDIMS 3
#define ARRAY_DIM (:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_4d_double_array_attribute
#define NDIMS 4
#define ARRAY_DIM (:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_5d_double_array_attribute
#define NDIMS 5
#define ARRAY_DIM (:,:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_6d_double_array_attribute
#define NDIMS 6
#define ARRAY_DIM (:,:,:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_7d_double_array_attribute
#define NDIMS 7
#define ARRAY_DIM (:,:,:,:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_scalar_double_attribute
#define NDIMS 0
#define SCALAR 1
#include "read_attribute.inc"
#undef SCALAR
#undef NDIMS
#undef SUB_NAME

#undef NATIVE_TYPE
#undef ARR_TYPE
#undef DATA_CLASS

!
! Routines to read n-dimensional string datasets
#define NATIVE_TYPE H5T_NATIVE_CHARACTER
#define ARR_TYPE CHARACTER(LEN=*)
#define DATA_CLASS H5T_STRING_F
#define STRING 1

#define SUB_NAME read_1d_string_array_attribute
#define NDIMS 1
#define ARRAY_DIM (:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_2d_string_array_attribute
#define NDIMS 2
#define ARRAY_DIM (:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_3d_string_array_attribute
#define NDIMS 3
#define ARRAY_DIM (:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_4d_string_array_attribute
#define NDIMS 4
#define ARRAY_DIM (:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_5d_string_array_attribute
#define NDIMS 5
#define ARRAY_DIM (:,:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_6d_string_array_attribute
#define NDIMS 6
#define ARRAY_DIM (:,:,:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_7d_string_array_attribute
#define NDIMS 7
#define ARRAY_DIM (:,:,:,:,:,:,:)
#include "read_attribute.inc"
#undef ARRAY_DIM
#undef NDIMS
#undef SUB_NAME

#define SUB_NAME read_scalar_string_attribute
#define NDIMS 0
#define SCALAR 1
#include "read_attribute.inc"
#undef SCALAR
#undef NDIMS
#undef SUB_NAME

#undef NATIVE_TYPE
#undef ARR_TYPE
#undef DATA_CLASS
#undef STRING
