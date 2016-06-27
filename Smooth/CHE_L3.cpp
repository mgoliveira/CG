/**
* @file    CHE_L3.cpp
* @author  Marcos Lage         <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner   <thomas.lewiner@polytechnique.org>
* @author  H�lio  Lopes         <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matm�dia
* @date    14/02/2006
*
* @brief  (Compact Half-Edge Structure - Level 3)
*/

#include <GL/glut.h>
// <GLUT/glut.h>
#include "CHE_L3.hpp"

using namespace std;

/*--------------------------------------------------------------------------*/
void CHE_L3::compute_CH()
/*--------------------------------------------------------------------------*/
/** Computes the vertex Boundary Curves table.*/
{
	cout << "CHE_L3::compute_CH..." ;

	vector<bool> vst;
	vst.resize(3*ntrig(), false);

	_CH.clear();
	for(HEid he=0; he<3*ntrig(); ++he)
	{
	   if(O(he) == -1 && vst[he] == false)
	   {
		   HEid he0 = he;
		   _CH.push_back(he0);  // add boundary curve representative
		   _ncurves++;          // increment number of curves
		   do                   // walk throght the boundary
		   {
				// _O[he0] = -(ncurves() + 1); 
				// set half--edge component
				// mark half--edge as visited
				vst[he0] = true; 
				while(O(next(he0)) >= 0) 
				{
					//get the next half--edge in the boundary
					he0 = O( next(he0) );
				}
				// real next boundary he.
				he0 = next(he0);
			}
			while(he0 != he);
		break;
	 }
	}
	vst.clear();
	cout << " done." << endl;
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
void CHE_L3::check()
/*--------------------------------------------------------------------------*/
/** Checks the mesh.*/
{
  CHE_L2::check();

  if( ncurves() != static_cast<int>(_CH.size()) ) {
    cout << "ncurves != _CH.size()." << endl; 
    return;
  }

  for(HEid i=0; i<3*ntrig(); ++i){
    if(O(i) < -(ncurves()+1)) {
      cout << "O(" << i << ") < ncurves()." << endl;
      return;
    }
  }

  for(int i=0; i<_ncurves; ++i)
	  if(O(_CH[i]) != -1) {
		  cout << "CH[" << i << "] = " << _CH[i] << " , Bound? = " << (O(_CH[i]) < 0) << endl;
		  return;
	  }
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
void CHE_L3::read_ply( const char* file )
/*--------------------------------------------------------------------------*/
/** Gets a boundary compound. Sets _CH*/
{
	
	
	
	// Stores Start && L3 time;
	clock_t start_time = static_cast<clock_t>(0.0),
			 L3_time = static_cast<clock_t>(0.0);

	start_time = clock();

	CHE_L2::read_ply( file );
	compute_CH();
	//CHE_L0::NLaplace(1);
	
	
	
	
	L3_time = clock();
	//cout << "L3 load time:" << static_cast<float>(L3_time-start_time)/static_cast<float>(CLOCKS_PER_SEC) << endl;
	cout << endl;
}
/*--------------------------------------------------------------------------*/
