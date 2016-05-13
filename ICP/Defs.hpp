#ifndef _DEFS_HPP_
#define _DEFS_HPP_


/** \brief Invalid integer*/


#include <map>
#include <set>
#include <queue>
#include <stack>
#include <ctime>
#include <cmath>
#include <vector>
#include <cfloat>
#include <iostream>
#include <algorithm>
#include <limits.h>
#include <armadillo>

#include "ply.h"

#define INV INT_MAX
#define DBG 0


using namespace std;

/** \brief Axis Definition */
enum{X,Y,Z};

/** Vertex id type; */
typedef int   Vid;

/** \brief Pet_triangle id*/
typedef int  TRid; 

/** \brief Pet_half-edge id*/
typedef int  HEid;

/** \brief Connected Compound id type */
typedef int   Cid;

typedef int OKid;
/*--------------------------------------------------------------------------*/
/** PLY Data*/
/*--------------------------------------------------------------------------*/
typedef struct PlyPoint
{
	float x, y, z;
	float nx, ny, nz;
}PlyPoint;
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/* list of property information for a PlyVertex */
static PlyProperty plyvert_props[]  = { 
  {"x",  Float32, Float32, offsetof( PlyPoint, x  ), 0, 0, 0, 0},
  {"y",  Float32, Float32, offsetof( PlyPoint, y  ), 0, 0, 0, 0},
  {"z",  Float32, Float32, offsetof( PlyPoint, z  ), 0, 0, 0, 0},
  {"nx", Float32, Float32, offsetof( PlyPoint, nx ), 0, 0, 0, 0},
  {"ny", Float32, Float32, offsetof( PlyPoint, ny ), 0, 0, 0, 0},
  {"nz", Float32, Float32, offsetof( PlyPoint, nz ), 0, 0, 0, 0},
};
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
typedef struct PlyFace
{
  unsigned char nverts;    /* number of Vertex indices in list */
  int *verts;              /* Vertex index list */
} PlyFace;
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/* list of property information for a PlyFace */
static PlyProperty plyface_props[]  = { 
{  "vertex_indices", Int32, Int32, offsetof( PlyFace,verts ),
   1, Uint8, Uint8, offsetof( PlyFace,nverts )  }
};
/*--------------------------------------------------------------------------*/
#endif
