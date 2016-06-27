/*
 * CHEL4.hpp
 *
 *  Created on: 9 de jun de 2016
 *      Author: marco
 */

#ifndef CHE_L4_HPP_
#define CHE_L4_HPP_

#include "CHE_L3.hpp"
#include <armadillo>

using namespace std;

class CHE_L4:public CHE_L3 {

protected:

	map<COTGid, float>   _COTG;
	map<ARid, float>   _Area;
	map<ARid, float>   _VoronoiArea;
	map<HEid, float>   _LT;
	map<LTid, float>   _LTV;

public:
	CHE_L4():CHE_L3(){};
	virtual ~CHE_L4(){ _V.clear(); _G.clear(); _O.clear(); _C.clear(); _VH.clear(); _EH.clear(); _CH.clear(); _COTG.clear(); _Area.clear(); _VoronoiArea.clear(); _LTV.clear();};

public:
	void compute_COTG();
	float compute_volume();
	void compute_Lenght();
	float compute_area(Vid);
	float voronoiArea(Vid);
	float edgeWeightCOTG(Vid , Vid);
	void scale_volume(float, float);
	vector<float> edgeWeightLGT(Vid);
	float vertexWeightCOTG(Vid );
	void laplacianSmooth();
	float* umbrella(Vid );
	float meancurvature();

public:
	/** \brief Checks the mesh*/
	   virtual void check();

public:
  /** \brief Reads a 3D model in the .ply format
    * \param file - const char* */
   virtual void  read_ply( const char* file );
};

#endif /* CHE_L4_HPP_ */
