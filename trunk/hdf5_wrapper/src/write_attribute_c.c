#include "hdf5.h"

void write_hdf5_attribute_(hid_t *pattr_id, hid_t *pmem_type_id, void *buf,
			 herr_t *err)
{
  *err = H5Awrite(*pattr_id,*pmem_type_id,buf);
}


