/*
 * CHEL4.cpp
 *
 *  Created on: 9 de jun de 2016
 *      Author: marco
 */

#include "CHE_L4.hpp"


/*--------------------------------------------------------------------------*/
void CHE_L4::compute_COTG()
/*--------------------------------------------------------------------------*/
/** Computes the cotangents in the triangle.*/
{
	cout << "CHE_L4::compute_COTG... " ;

	_COTG.clear();

	for(int v = 0; v < nvert(); ++v){

		//For each vertex gets the first incident half--edge

		HEid h = VH(v), h1 = h, h0 = h;
		Vid ia, ib, ic;

		//Walks through vertices in the star and compute the cotangent.

		ia = v;

		if(v_valid(ia)){

			do{
				h1 = h;
				h  = O(prev(h1));

				TRid t =  trig(h1);

				ib = V(next(h1));
				ic = V(prev(h1));

				float Px =  G(ib).x() - G(ia).x(); float Py = G(ib).y() - G(ia).y(); float Pz = G(ib).z() - G(ia).z();
				float Qx =  G(ic).x() - G(ia).x(); float Qy = G(ic).y() - G(ia).y(); float Qz = G(ic).z() - G(ia).z();

				float P_dot_Q = Px * Qx + Py * Qy + Pz * Qz;

				float Norm_P = sqrt(Px * Px + Py * Py + Pz * Pz);
				float Norm_Q = sqrt(Qx * Qx + Qy * Qy + Qz * Qz);

				double k = P_dot_Q / (Norm_P * Norm_Q);

				// k should belong to interval (-0.99, 0.99)

				float theta = acos( max( -0.99, min( 0.99, k ) ) );

				_COTG.insert((pair<COTGid, float>(make_pair(ia, t), 1.0/tan(theta))));

			}
			while( (h != -1) && (h != h0) );
		}

	}

	cout << "done." << endl;

}

/*-------------------------------------------------------------------------- */
float CHE_L4::vertexWeightCOTG(Vid v)
/*-------------------------------------------------------------------------- */
/* Compute the cotangent weight in vertex.*/
{
	float weight = 0.0;
	vector<Vid> vertex = R_00(v);

	for(vector<Vid>::iterator it = vertex.begin(); it != vertex.end(); ++it){
		weight += edgeWeightCOTG(v, *it);
	}

	return weight;
}

/*--------*/

/*--------------------------------------------------------------------------*/
float CHE_L4::edgeWeightCOTG(Vid vo, Vid v)
/*--------------------------------------------------------------------------*/
/** Compute the weight in edge .*/
{
	vector<Vid> starVertex = R_00(vo);

	double weight = 0.0;

	HEid hf = -1;

	if ( find(starVertex.begin(), starVertex.end(), v) != starVertex.end() ){

		for (map<HEid, Edge>::iterator it = _EH.begin(); it != _EH.end();  ++it)
		{
			if ( (it->second.a() == vo && it->second.b() == v)  || (it->second.a() == v && it->second.b() == vo) )
			{
				hf = it->first;
				break;
			}
		}

		if((O(hf) != -1)){
			COTGid key_01 = make_pair(V(prev(O(hf))), trig(O(hf)));
			weight += _COTG[key_01];
		}

		COTGid key_02 = make_pair(V(prev(hf)), trig(hf));
		weight += _COTG[key_02];

	}

	return max(weight/2.0, 0.0);

}

/*--------------------------------------------------------------------------*/
float CHE_L4::voronoiArea(Vid v)
/*--------------------------------------------------------------------------*/
/** Compute the voronoi area of star of the vertex of mesh.*/
{

	float area = 0.0;

	//Gets the first incident half--edge

	HEid h = VH(v), h1, h0 = h;
	Vid ia, ib, ic;

	//Walks through vertices in the triangle and compute the voronoi area of triangle.

	ia = v;

	if(v_valid(ia)){

		do{
			h1 = h;

			h  = O(prev(h1));

			ib = V(next(h1));
			ic = V(prev(h1));

			float t = trig(h1);

			float Px =  G(ib).x() - G(ia).x(); float Py = G(ib).y() - G(ia).y(); float Pz = G(ib).z() - G(ia).z();
			float Qx =  G(ic).x() - G(ia).x(); float Qy = G(ic).y() - G(ia).y(); float Qz = G(ic).z() - G(ia).z();

			Vertex Cross_product;

			Cross_product.set_x( (Py * Qz) - (Qy * Pz) );
			Cross_product.set_y( (Pz * Qx) - (Qz * Px) );
			Cross_product.set_z( (Px * Qy) - (Qx * Py) );

			float a = (1.0/2.0) * Cross_product.norm(Cross_product.x(), Cross_product.y(), Cross_product.z());

			float norm_P = Px * Px + Py * Py + Pz * Pz;
			float norm_Q = Qx * Qx + Qy * Qy + Qz * Qz;

			COTGid key_01 = make_pair(V(h1), trig(h1));
			COTGid key_02 = make_pair(V(prev(h1)), trig(h1));
			COTGid key_03 = make_pair(V(next(h1)), trig(h1));

			float A = _COTG[key_01];
			float B = _COTG[key_02];
			float C = _COTG[key_03];

			if( A < 0 || B < 0 || C <0 ){
				if(A < 0){
					area += 0.5* a ;
				}else area += 0.25* a ;
			}else area += (1.0/8.0) * ( norm_P * B + norm_Q * C );

		}while( (h != -1) && (h != h0) );

	}

	return area;
}


/*--------------------------------------------------------------------------*/
void CHE_L4::scale_volume(float Voli, float Voln)
/*--------------------------------------------------------------------------*/
/** Aplly the scale volume of mesh.*/
{

	float beta = pow(Voli/Voln, 1.0/3.0);
	for(int i = 0; i < (int) _G.size(); ++i){

		G(i).set_x( beta * G(i).x() );
		G(i).set_y( beta * G(i).y() );
		G(i).set_y( beta * G(i).y() );

	}

}


/*--------*/

/*--------------------------------------------------------------------------*/
float CHE_L4::compute_volume()
/*--------------------------------------------------------------------------*/
/** Compute the volume of mesh.*/
{

	float volume = 0.0;

	for(int v = 0; v < nvert(); ++v){

		//For each vertex gets the first incident half--edge

		HEid h = VH(v);
		Vid ia, ib, ic;

		//Walks through vertices in the triangle and compute the volume of triangle.

		ia = v;

		if(v_valid(ia)){

			ib = V(next(h));
			ic = V(prev(h));

			float Px =  G(ib).x() - G(ia).x(); float Py = G(ib).y() - G(ia).y(); float Pz = G(ib).z() - G(ia).z();
			float Qx =  G(ic).x() - G(ia).x(); float Qy = G(ic).y() - G(ia).y(); float Qz = G(ic).z() - G(ia).z();

			Vertex Cross_product;

			Cross_product.set_x( (Py * Qz) - (Qy * Pz) );
			Cross_product.set_y( (Pz * Qx) - (Qz * Px) );
			Cross_product.set_z( (Px * Qy) - (Qx * Py) );

			Vertex g;

			g.set_x( (1.0/3.0) * ( G(ia).x() + G(ib).x() + G(ic).x() ));
			g.set_y( (1.0/3.0) * ( G(ia).y() + G(ib).y() + G(ic).y() ));
			g.set_z( (1.0/3.0) * ( G(ia).z() + G(ib).z() + G(ic).z() ));

			volume +=  Cross_product.x() * g.x() + Cross_product.y() * g.y() + Cross_product.z() * g.z();

		}

	}

	return volume/6.0;

}

/*--------------------------------------------------------------------------*/
void CHE_L4::laplacianSmooth(float Voli, float l)
/*--------------------------------------------------------------------------*/
/** Laplacian smooth.*/
{

	arma::mat Bx(nvert(), 1), By(nvert(), 1), Bz(nvert(), 1);

	arma::mat I(nvert(), nvert()); I.eye();

	arma::mat D(nvert(), nvert()); D.zeros();

	arma::mat A(nvert(), nvert()); A.zeros();

	arma::mat L;

	for(int i = 0; i < nvert(); ++i){

		Bx(i, 0) = G(i).x();
		By(i, 0) = G(i).y();
		Bz(i, 0) = G(i).z();

		for(int j = 0; j < nvert(); ++j){

			if( i != j ){

				A(i, j) = edgeWeightCOTG(i, j);

			}else{

				D(i, j) = 1.0 / (voronoiArea(i));
				A(i, j) = -vertexWeightCOTG(i);

			}

		}

	}

	float min[3], max[3];
	bounding_box(min, max);

	float x = max[0] -min[0];
	float y = max[1] -min[1];
	float z = max[2] -min[2];

	float factor_norm = sqrt(x * x + y * y + z * z) * 0.01;

	float dt = l * factor_norm;

	L =  D * A;

	arma::sp_mat K;

	K = I - dt*L;

	arma::mat X = arma::spsolve( K, Bx );
	arma::mat Y = arma::spsolve( K, By );
	arma::mat Z = arma::spsolve( K, Bz );

	for(int i = 0; i < nvert(); ++i){

		G(i).set_x(X(i,0));
		G(i).set_y(Y(i,0));
		G(i).set_z(Z(i,0));

	}

	float Volf = compute_volume();

	scale_volume(Voli, Volf);

	cout << "Smooth... done" << endl;

}

/*--------*/

/*--------------------------------------------------------------------------*/
void CHE_L4::check()
/*--------------------------------------------------------------------------*/
/** Checks the mesh.*/
{

	CHE_L3::check();

	if(3*ntrig() != static_cast<int>(_COTG.size()) ){
		cout << "CHE_L4:: Error 3*ntrig()!= COTG.size" << endl;
		return;
	}

	for(map<COTGid, float>::iterator it = _COTG.begin(); it != _COTG.end(); ++it){

		Vid v = (*it).first.first;

		HEid h = VH(v);

		TRid t =  (*it).first.second;

		if(!v_valid(v)){
			cout << "CHE_L4:: Error vertex " << v << " invalid" << endl;
			return;
		}

		if(!he_valid(h)){
			cout << "CHE_L4:: Error half-edge " << h << " invalid" << endl;
			return;
		}

		if(!tr_valid(t)){
			cout << "CHE_L4:: Error triangle " << (*it).first.second << " invalid" << endl;
			return;
		}

		if(isnanf((*it).second)){
			cout << "CHE_L4:: Error COTG. Angle NaN in vertex " << (*it).first.first << "and in triangle " << t << endl;
			return;
		}

	}
}
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
void CHE_L4::read_ply( const char* file )
/*--------------------------------------------------------------------------*/
/** Gets a boundary compound. */
{

	// Stores Start && L3 time;
	clock_t start_time = static_cast<clock_t>(0.0), L4_time = static_cast<clock_t>(0.0);

	start_time = clock();

	CHE_L3::read_ply( file );
	compute_COTG();

	L4_time = clock();
	//cout << "L3 load time:" << static_cast<float>(L3_time-start_time)/static_cast<float>(CLOCKS_PER_SEC) << endl;
	cout << endl;
}
/*--------------------------------------------------------------------------*/


