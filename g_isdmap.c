/* 
 * 
 * Tim Connolly - tconnolly@ucmerced.edu
 * Copyright (c) 2014, Regents of the University of California
 * Released under BSD 2-Clause License (see "LICENSE" file)
 * 
 * This code was modified from the file src/tools/g_gyrate.c
 */


/* This is just a wrapper binary.
* The code that used to be in g_shape.c is now in gmx_shape.c,
* where the old main function is called gmx_shape().
*/
int
main(int argc, char *argv[])
{
  gmx_isdmap(argc,argv);
  return 0;
}


  
