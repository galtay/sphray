
/* Return 0 for big endian, 1 for little endian */
void byte_order_(int *i)
{                                                   
  int   one = 1;
  char* endptr = (char *) &one;
  *i = (*endptr);
  return;
} 
