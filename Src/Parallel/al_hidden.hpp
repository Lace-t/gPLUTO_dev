////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Internal include file for the ArrayLib.

  Internal include file for the ArrayLib.
  If this is defined, each node has a copy of the indexes from all 
  other nodes. This is mainly useful to deal with non-uniform array 
  distributions.

  \author A. Malagoli (University of Chicago)
  \author G. Muscianisi (g.muscianisi@cineca.it)

  \date   Aug 29, 2012
*/
////////////////////////////////////////////////////////////////////////
#ifndef __AL_HIDDEN
#define __AL_HIDDEN


#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "al_defs.hpp"

/* Include AL macros and definitions */
#include "al_codes.hpp"

/* Include functions prototypes */
#include "al_proto.hpp" 


/* Type definition for the SZ structure */
typedef struct szz{
  int compiled;         /* AL_TRUE if this is a compiled array */
  MPI_Datatype type;    /* The elementary data type of the array
                               e.g.: MPI_FLOAT, MPI_INT .. */
  int type_size;        /* Size of the elmentary data type */
  int ndim;             /* Number of array dimensions (max: 5) */
  int rank;             /* Rank of this portion of the array */
  int size;             /* Size of the communicator which the 
                               array is associated with */
  MPI_Comm comm;        /* The communicator which the array is
			   associated with */
  MPI_Comm cart_comm;   /* The derived cartesian communicator */
  MPI_Comm oned_comm[AL_MAX_DIM]; /* 1D communicators along directions */
      /* (Useful for performing global operations along single directions */
  long int buffsize;    /* Size of the local buffer for the array */
  int arrdim[AL_MAX_DIM];  /* Global dimensions of the array */
  int larrdim[AL_MAX_DIM]; /* Local dimensions of the array
                              in this dimension WITHOUT ghost points */
  int larrdim_gp[AL_MAX_DIM];  /* Local size of array in this direction
                                  WITH ghost points (compare to larrdim)*/
  int beg[AL_MAX_DIM];     /* Global address of beginning of array
			      in this dimension (C-convention) */
  int end[AL_MAX_DIM];     /* Global address of ending of array
			      in this dimension (C-convention) */

  int *begs, *ends;        /* These pointers  can be used to 
			      gather the information about the 
                              distribution of the indexes across 
                              all the nodes */

  int lbeg[AL_MAX_DIM];    /* Local address of beginning of array
			      in this dimension (C-convention) */
  int lend[AL_MAX_DIM];    /* Local address of ending of array
			      in this dimension (C-convention) */
  int bg[AL_MAX_DIM];      /* Size of ghost point at beginning of 
			      array in this dimension */
  int eg[AL_MAX_DIM];      /* Size of ghost point at beginning of 
			      array in this dimension */
  int offset[AL_MAX_DIM];  /* Local array offsets for local do loops */
  int stride[AL_MAX_DIM];  /* We will replace offset with stride EVENTUALLY */
  int isparallel[AL_MAX_DIM];  /* AL_TRUE if the array is distributed 
                                  along this dimension [Default: AL_TRUE]*/
  int isperiodic[AL_MAX_DIM];  /* AL_TRUE if the array is periodic in 
                                  this dimension [Default: AL_TRUE] */
  int isstaggered[AL_MAX_DIM]; /* AL_TRUE if the array is staggered in 
                                  this dimension [Default: AL_FALSE]*/
  int left[AL_MAX_DIM];  /* Rank of left node in this dimension in the
                            cartesian communicator topology */
  int right[AL_MAX_DIM]; /* Rank of right node in this dimension in the
                            cartesian communicator topology */
  int lrank[AL_MAX_DIM]; /* Rank of this node along this dimension of 
                            the cartesian nodes topology */
  int lsize[AL_MAX_DIM]; /* Number of nodes along this dimension of
                            the cartesian nodes topology */
  int sendb1[AL_MAX_DIM]; /* Send buffer pointer for the SendRecv operation */
  int sendb2[AL_MAX_DIM];
  int recvb1[AL_MAX_DIM]; /* Receive buffer pointers for the SendRecv operation */
  int recvb2[AL_MAX_DIM];
  int tag1[AL_MAX_DIM], tag2[AL_MAX_DIM]; /* Tags for SendRecv operation */
  MPI_Datatype  strided[AL_MAX_DIM]; /* Strided data type corresponding to
                                        the ghost point buffer in this 
                                        dimension */
  /* We add two separate strided data types in order to allow for 
     a different exchange of ghost points from Left to Right than
     from Right to Left. This is needed when we exchange boundary
     points in strided or overlapping data types */
  MPI_Datatype  type_lr[AL_MAX_DIM]; /* Strided data type for the 
                                        Left->Right exchange */
  MPI_Datatype  type_rl[AL_MAX_DIM]; /* Strided data type for the
                                        Right->Left exchange */
  int pio_offset;       /* The offset for PIO */
  
  MPI_Datatype gsubarr; /* Global subarray for MPI-IO */
  MPI_Datatype lsubarr; /* Local subarray for MPI-IO */

  MPI_Datatype gsubarr_stag[AL_MAX_DIM]; /* Global staggered subarrays for MPI-IO */
  MPI_Datatype lsubarr_stag[AL_MAX_DIM]; /* Local staggered subarrays for MPI-IO */

  MPI_Offset io_offset;  /* Offset used to store file pointer */
  MPI_File ifp;          /* Pointer to the file this array is to be
                            written using MPI-IO */
} SZ;

#endif /* End ifndef __AL_HIDDEN */
