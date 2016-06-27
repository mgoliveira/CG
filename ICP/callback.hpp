/*
 * callback.hpp
 *
 *  Created on: 5 de abr de 2016
 *      Author: marco
 */

#ifndef CALLBACK_HPP_
#define CALLBACK_HPP_

#include<vector>

const int k = 3;

struct Node{
    float point[6];
    Node *left, *right;
};

struct NewData{
	std::vector<float*> Malha;
	float erro;
};

#include "CHE_L0.hpp"
#include <armadillo>
void reshape(int , int );
void motion(int , int );
void mouse(int , int , int , int );
void keyboard(unsigned char , int , int );
void zoonIn(void);
void zoonOut(void);
void trackball(int , int , int , int , float* );
void display(void);
void init(void);
void drawDot(float, float, float);
void reset(void);
void axe(void);
void plano(void);
float* XwinToView(int, int);
void curvaBS(void);
void drawLine(int);
void grid(void);

vector<float*> vectorDataICP(CHE_L0, int);
arma::mat covariance(vector<float*>, vector<float*>, float*, float*);
arma::mat rotation(arma::mat );
arma::mat translation(arma::mat ,  Vertex , Vertex );
arma::mat reflection(arma::mat );
float* nearpoint(float* , Node* );
float* VertexToFloat(CHE_L0 , int );
float* centroid(vector<float*>);
struct Node* kdtree(vector<float*>);
struct Node* CHEkdtree(CHE_L0);
void ICPPointPoint(CHE_L0, CHE_L0, float);
float* ArmtoFloat(arma::mat , int , int );
float* ArmtoRot(arma::mat);
vector<float*> CHEtoVector(CHE_L0 );
float dis( float* ,  float* );
void viewMesh();
void multiView();
void boudingSphere();

void TransCent(CHE_L0 , CHE_L0 );

arma::mat PointToPlane(vector<float*> , vector<float*> );

void ICPPointPlane(CHE_L0, CHE_L0, float);


Node* insertR(Node*, vector<float*>& , unsigned );

#endif /* CALLBACK_HPP_ */
