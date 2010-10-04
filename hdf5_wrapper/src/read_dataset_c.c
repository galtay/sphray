#include "hdf5.h"
// Routine to read an HDF5 dataset from Fortran - this should allow reading
// of data types not defined in the HDF5 Fortran interface, although it does
// circumvent some Fortran error checking...

void read_hdf5_dataset_(hid_t *pdataset_id,hid_t *pmem_type_id,
			hid_t *pmemspace_id, hid_t *pdspace_id,
			void * buf,herr_t *err)
{
  *err = H5Dread( *pdataset_id, *pmem_type_id, *pmemspace_id,
		  *pdspace_id, H5P_DEFAULT, buf);
}


