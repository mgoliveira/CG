/*
 * callback.cpp
 *
 *  Created on: 5 de abr de 2016
 *      Author: marco
 */

#include "CHE_L0.hpp"

#define bool int
#define false 0
#define true 1
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <math.h>
#include <iostream>
#include "callback.hpp"
#include <armadillo>
#include <gsl/gsl_math.h>

float 	dz = 45,
		dx = 0,
		dy = 0,
		theta = 0,
		phi=0,
		raio = 45,
		best_dis = FLT_MAX,
		scl = 10,
		lastPos[3] = {0.0, 0.0, 0.0};

bool isEnableTrackBall = false,
	 isViewPrincipal = true,
	 isBS = true;

GLfloat Transform[16] = {1.0f,  0.0f,  0.0f,  0.0f,
		                 0.0f,  1.0f,  0.0f,  0.0f,
		                 0.0f,  0.0f,  1.0f,  0.0f,
		                 0.0f,  0.0f,  0.0f,  1.0f
};

CHE_L0 objeto1;
CHE_L0 objeto2;

void init(){

	glClearColor(0.3, 0.3, 0.3, 1.0);
	glPointSize(5.0);
	glViewport(0, 0, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(dz, 1.5, 0.01, 300);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glEnable(GL_DEPTH_TEST);

	objeto1.read_ply("./PlyFiles/b06_low.ply"); // Lê a malha a ser ajustada a referência
	objeto2.read_ply("./PlyFiles/b07_low.ply"); // Lê a malha referência

	ICP(objeto1, objeto2, 0.001); // Chama a função ICP. O último parâmetro é a torelância do ajuste
}

void display(void){

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	dx = raio*cos(phi)*sin(theta);
	dy = raio*sin(phi)*sin(theta);
	dz = raio*cos(theta);

	if(isViewPrincipal == true){
		glLoadIdentity();
		gluLookAt(dx, dy, dz, 0.0, 0.0, 0.0, 0,1,0);
		viewMesh();
	}
	else{
		multiView();
	}

	glutSwapBuffers();

}

void multiView(){

	glViewport(0, 0, glutGet(GLUT_WINDOW_WIDTH)/2, glutGet(GLUT_WINDOW_HEIGHT)/2);
	glLoadIdentity();
	gluLookAt(0.0, 0.0, -45, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	viewMesh();

	glViewport(glutGet(GLUT_WINDOW_WIDTH)/2, 0, glutGet(GLUT_WINDOW_WIDTH)/2, glutGet(GLUT_WINDOW_HEIGHT)/2);
	glLoadIdentity();
	glColor3f(0.1, 0.2, 0.2);
	gluLookAt(0.0, 0.0, 45, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	viewMesh();
	glColor3f(0.0, 0.0, 0.0);

	glViewport(0, glutGet(GLUT_WINDOW_HEIGHT)/2, glutGet(GLUT_WINDOW_WIDTH)/2, glutGet(GLUT_WINDOW_HEIGHT)/2);
	glLoadIdentity();
	gluLookAt(0.0, 45, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
	viewMesh();

	glViewport(glutGet(GLUT_WINDOW_WIDTH)/2, glutGet(GLUT_WINDOW_HEIGHT)/2, glutGet(GLUT_WINDOW_WIDTH)/2, glutGet(GLUT_WINDOW_HEIGHT)/2);
	glLoadIdentity();
	gluLookAt(0.0, -45, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
	viewMesh();
}

void viewMesh(){

		if(isViewPrincipal == true) {
			glColor3f(0.0, 0.0, 0.0);
		}

		grid();

		axe();

		glMultMatrixf( (GLfloat *) Transform);

		glColor3f(0.0, 0.0, 0.0);

		glColor3f(1.0, 0.0, 0.0);
		objeto1.draw_wire();

		glColor3f(0.0, 1.0, 0.0);
		objeto2.draw_wire();

		if(isBS){
			glPushMatrix();
			glScalef(scl, scl, scl);

			boudingSphere();
			glPopMatrix();
		}
}
void boudingSphere(){

	float k = 0;

	float h = 0;

	float r = 1;

	float x = 0, y = 0, z = 0;

	glLineWidth(1);

	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	for (int i = 0; i < 180; i++)
	{
		x = r * cos(i) - h;
		y = r * sin(i) + k;
		glVertex3f(x + k, y - h,0);

		x = r * cos(i + 0.1) - h;
		y = r * sin(i + 0.1) + k;
		glVertex3f(x + k,y - h,0);
	}
	glEnd();

	y = 0; 	z = 0;
	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	for (int i = 0; i < 180; i++)
	{
		y = r * cos(i) - h;
		z = r * sin(i) + k;
		glVertex3f(0,y + k, z - h);

		y = r * cos(i + 0.1) - h;
		z = r * sin(i + 0.1) + k;
		glVertex3f(0, y + k, z - h);
	}
	glEnd();

	x = 0; 	z = 0;
	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINES);
	for (int i = 0; i < 180; i++)
	{
		x = r * cos(i) - h;
		z = r * sin(i) + k;
		glVertex3f(x + k, 0, z - h);

		x = r * cos(i + 0.1) - h;
		z = r * sin(i + 0.1) + k;
		glVertex3f(x + k, 0, z - h);
	}
	glEnd();

	glLineWidth(1);
	glColor3f(0.0, 0.0, 0.0);

}

void reshape(int w, int h){
	if (h == 0) h = 1;

	glMatrixMode (GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, (GLsizei) w, (GLsizei) h);
	gluPerspective(dz,(GLfloat) w/(GLfloat) h, 0.01, 300);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void mouse(int button, int state, int x, int y){

	if(button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN){
		isEnableTrackBall = true;
		trackball(x, y, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT), lastPos);
	}else isEnableTrackBall = false;

	if ((button == 3) || (button == 4)){
		if (state == GLUT_UP) return;
		if(button == 3){
			zoonIn();
		}
		else zoonOut();
	}

}

void motion(int x, int y){

	float* curPos = new float[3]();
	float* axis   = new float[3]();

	if(isEnableTrackBall){

		trackball(x, y, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT), curPos);

		float dot = curPos[0]*lastPos[0] + curPos[1]*lastPos[1] + curPos[2]*lastPos[2];

		float angle = acos(MIN(1.0f, dot));

		axis[0] = lastPos[1]*curPos[2] - lastPos[2]*curPos[1];
		axis[1] = lastPos[2]*curPos[0] - lastPos[0]*curPos[2];
		axis[2] = lastPos[0]*curPos[1] - lastPos[1]*curPos[0];

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glRotatef(angle*180/M_PI, axis[0], axis[1], axis[2]);
		glMultMatrixf( (GLfloat *) Transform );
		glGetFloatv(GL_MODELVIEW_MATRIX,(GLfloat *) Transform);

		lastPos[0] = curPos[0];
		lastPos[1] = curPos[1];
		lastPos[2] = curPos[2];
	}

	delete(curPos);
	delete(axis);

	glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y){
	switch (key) {
	case 'a':
		theta = theta + 0.1;
		break;
	case 'b':
		if(isBS){
			isBS = false;
		}else
			isBS = true;
		break;
	case 'p':
		if(isViewPrincipal){
			isViewPrincipal = false;
		}else
			isViewPrincipal = true;
			reshape(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
		break;
	case 'd':
		theta = theta - 0.1;
		break;
	case 'w':
		phi = phi + 0.1;
		break;
	case 's':
		phi = phi - 0.1;
		break;
	case '+':
		scl = scl + 0.1;
		break;
	case '-':
		scl = scl - 0.1;
		break;
	}
	glutPostRedisplay();
}


void zoonIn(void){
	if(dz < 180){
		raio += 0.5;
		glutPostRedisplay();
	}
}

void zoonOut(void){
	if(dz > 0){
		raio -=0.5;
		glutPostRedisplay();
	}
}

void trackball(int x, int y, int width, int height, float* v)
{
	v[0] = (2.0*x / width - 1.0);
	v[1] = -1*(2.0 * y / height - 1.0);
	v[2] = 0;
	float d = sqrt(v[0]*v[0] + v[1]*v[1]);
	float a = 1.0 / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	if (d <= 1.0){
		v[2] = sqrt(1.0 - d);
		a = 1.0 / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
		v[0] *= a; v[1] *= a; v[2] *= a;
	}
	else{
		v[0] *= a; v[1] *= a; v[2] *= a;
	}
}

void axe(){

	glPushMatrix();

	glTranslatef(15,-6,15);
	glScalef(0.3, 0.3, 0.3);
	glMultMatrixf( (GLfloat *) Transform);

	glColor3f(1.0, 0.0, 0.0);

	glBegin (GL_LINES);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(5.0, 0.0, 0.0);
	glEnd();

	glPushMatrix();
	glTranslated(5.0,0.0, 0.0);
	glRotated(90, 0.0, 1.0, 0.0);
	glutWireCone(0.3, 0.7, 15, 15);
	glPopMatrix();

	glColor3f(0.0, 1.0, 0.0);
	glBegin (GL_LINES);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 5.0, 0.0);
	glEnd();

	glPushMatrix();
	glTranslated(0.0, 5.0, 0.0);
	glRotated(90, -1.0, 0.0, 0.0);
	glutWireCone(0.3, 0.7, 15, 15);
	glPopMatrix();

	glColor3f(0.0, 0.0, 1.0);
	glBegin (GL_LINES);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, 5.0);
	glEnd();

	glPushMatrix();
	glTranslated(0.0, 0.0, 5.0);
	glRotated(90, 0.0, 0.0, 1.0);
	glutWireCone(0.3, 0.7, 15, 15);
	glPopMatrix();

	glPopMatrix();

}

void grid(){

	glBegin(GL_LINES);
	for (GLfloat i = -20.0; i <= 20.0; i += 1) {
		glVertex3f(i, -10, 20.0); glVertex3f(i, -10, -20.0);
		glVertex3f(20.0, -10, i); glVertex3f(-20.0, -10, i);
	}

	glEnd();
}

//Cria uma estrutura para armazenar os pontos selecionados de uma malha para realizar o ICP

vector<float*> vectorDataICP(CHE_L0 malha, int num){

	vector<float*> Mesh = vector<float*>() ;


	for(int i = 0; i < num; ++i){
		Mesh.push_back(VertexToFloat(malha, i));

	}
	return Mesh;
}

//Calcula o centróide de um conjunto de dados

float* centroid( vector<float*> data){

	float* cent = new float[3]();

	cent[0] = 0;
	cent[1] = 0;
	cent[2] = 0;

	for(vector<float*>::iterator it = data.begin(); it != data.end(); ++it)
		{
		cent[0] += (*it)[0];
		cent[1] += (*it)[1];
		cent[2] += (*it)[2];
		}

	cent[0] = cent[0]/(float)data.size();
	cent[1] = cent[1]/(float)data.size();
	cent[2] = cent[2]/(float)data.size();

	return cent;

}

//Calcula a matriz de covariância

arma::mat covariance(vector<float*> datA, vector<float*> datB, float* centA,  float* centB){

	arma::mat A(3,1), B(3,1), H(3,3), ca(3,1), cb(3,1);

	H.zeros();
	A.zeros();
	B.zeros();

	vector<float*>::iterator itB = datB.begin();

	for(vector<float*>::iterator itA = datA.begin(); itA != datA.end(); ++itA)
	{
		A(0,0) = (float) (*itA)[0] - centA[0];
		A(1,0) = (float) (*itA)[1] - centA[1];
		A(2,0) = (float) (*itA)[2] - centA[2];

		B(0,0) = (float) (*itB)[0] - centB[0];
		B(1,0) = (float) (*itB)[1] - centB[1];
		B(2,0) = (float) (*itB)[2] - centB[2];

		H += A*B.t();

		++itB;

	}

	return H;
}

arma::mat rotation(arma::mat G){

	arma::mat U, V;
	arma::vec s;

	U.zeros();
	V.zeros();
	s.zeros();

	arma::svd(U,s,V,G);

	arma::mat H = V*U.t();

	return V*U.t();
}

// Calcula a matriz de translação

arma::mat translation( arma::mat R,  float* centA,  float* centB){

	arma::mat CentA(3,1), CentB(3,1);

	CentA.zeros();
	CentB.zeros();

	CentA(0,0) = centA[0];
	CentA(1,0) = centA[1];
	CentA(2,0) = centA[2];

	CentB(0,0) = centB[0];
	CentB(1,0) = centB[1];
	CentB(2,0) = centB[2];

	return (-1.0)*R*CentA + CentB;
}

//Corrige a reflexão se houver

arma::mat reflection(arma::mat R){

	if(arma::det(R) < 0){
		for(int i = 0; i < (int) R.n_rows; i++){
			R(i,2) = R(i,2)*(-1);
		}
	}

	return R;
}

//Converte um vertex de uma malha em um array tipo float

float* VertexToFloat(CHE_L0 malha, int i){

	float* ponto = new float[6]();

	ponto[0] = malha.G(i).x();
	ponto[1] = malha.G(i).y();
	ponto[2] = malha.G(i).z();
	ponto[3] = malha.G(i).nx();
	ponto[4] = malha.G(i).ny();
	ponto[5] = malha.G(i).nz();

	return ponto;
}

//Constrói uma KdTree a partir de uma malha

struct Node* CHEkdtree(CHE_L0 data){

	struct Node *root = NULL;

	vector<float*> dados = vector<float*>();

	for(int i = 0; i < data.nvert(); i++){
		dados.push_back(VertexToFloat(data, i));
	}

	root = insertR(root, dados, 0);

	return root;
}

//Constrói um novo Nó

struct Node* newNod(float* arr)
{
    struct Node* temp = new Node;

    for (int i=0; i<k; i++)
        temp->point[i] = arr[i];

    temp->left = NULL;
    temp->right = NULL;

    return temp;
}

//Determina onde será inserido na KdTree o novo ponto

Node* insertR(Node* root, vector<float*> & dados, unsigned depth)
{
	unsigned cd = depth % k;

	struct compare{
		int cd;
		compare(int cmd){this->cd = cmd;}
		bool operator()(float* a, float* b) {return (a[cd] < b[cd]);}
	};

	size_t half_size = dados.size() / 2;

	sort(dados.begin(), dados.end(), compare(cd));

	float* median = *(dados.begin() + half_size);

	if(root == NULL){
		root = newNod(median);
	}

	vector<float*> rootLeft =  vector<float*>();
	vector<float*> rootRight =  vector<float*>();

	if(half_size > 0){
		for(vector<float*>::iterator it = dados.begin(); it != dados.begin() + half_size; ++it){
			rootLeft.push_back(*it);
		}

		for(vector<float*>::iterator it = dados.begin() + half_size; it != dados.end(); ++it){
			rootRight.push_back(*it);
		}
		root->left = insertR(root->left, rootLeft, depth + 1);
		root->right = insertR(root->right, rootRight, depth + 1);
	}

	return root;
}

//Calcula a distância entre dois pontos

float dis( float* point1,  float* point2){
	float dist = 0.0;
	for(int i = 0; i < k; ++i){
		dist += (point1[i]-point2[i])*(point1[i]-point2[i]);
	}
	return dist;
}

//Determina o ponto mais próximo na KdTree

float* neareast(Node *root, float* point, int depth) {

	int cd = depth % k;

	float dx = point[cd] - (root->point)[cd];

	Node *near = dx < 0 ? root->left  : root->right;
	Node *far  = dx < 0 ? root->right : root->left;

	float* best = (near == NULL) ? root->point : neareast(near, point, depth + 1);

	if(dis(root->point, point) < dis(best, point)){
		best = root->point;
	}

	if(far != NULL){
		if(pow((root->point[cd] - point[cd]),2) < dis(best,point)){
			float* otherBest = neareast(far, point, depth + 1);
			if(dis(otherBest, point) < dis(best, point)){
				best = otherBest;
			}
		}

	}
	return best;
}

//Determina o ponto mais próximo de um dado ponto

float* nearpoint(float* ponto, Node* KdTree){

	return neareast(KdTree, ponto,0);
}

// Função ICP

void ICP(CHE_L0 malha1, CHE_L0 malha2, float threshold){

	struct Node* rootA = NULL;

	rootA = CHEkdtree(malha2);

	vector<float*> Mesh = vector<float*>();

	float error = 0.0, preverror = 0.0;

	bool target = false;

	int i = 0;

	Mesh = vectorDataICP(malha1, malha1.nvert());

	while(!target){

		vector<float*> MeshAux = vector<float*>();
		vector<float*> matchMesh = vector<float*>();

		for(vector<float*>::iterator it = Mesh.begin(); it != Mesh.end(); ++it){

			float* ponto = nearpoint(*it, rootA);

			float d = dis(*it, ponto);

			matchMesh.push_back(ponto);

			preverror += d;
		}

		preverror = preverror/(float) Mesh.size();

		float* centA = centroid(Mesh);

		float* centB = centroid(matchMesh);

		arma::mat H = covariance(Mesh, matchMesh, centA, centB);

		arma::mat R = rotation(H);

		R = reflection(R);

		arma::mat T = translation(R, centA, centB);

		//Realização dos ajustes e determinação do novo erro

		arma::mat newpoint(3,1);
		arma::mat newnorma(3,1);
		int j = 0;
		newpoint.zeros(); newnorma.zeros();

		for(vector<float*>::iterator it = Mesh.begin(); it != Mesh.end(); ++it){

			float* newponto = new float[6]();

			newpoint(0,0) = (*it)[0] ; newpoint(1,0) = (*it)[1]; newpoint(2,0) = (*it)[2];
			newnorma(0,0) = (*it)[3] ; newnorma(1,0) = (*it)[4]; newnorma(2,0) = (*it)[5];

			newpoint = R*newpoint;
			newnorma = R*newnorma;
			newpoint = newpoint + T;
			newnorma = newnorma + T;

			objeto1.G(j).set_x(newpoint(0,0));objeto1.G(j).set_y(newpoint(1,0));objeto1.G(j).set_z(newpoint(2,0));
			objeto1.G(j).set_nx(newnorma(0,0));objeto1.G(j).set_ny(newnorma(1,0));objeto1.G(j).set_nz(newnorma(2,0));

			newponto[0] = newpoint(0,0);newponto[1] = newpoint(1,0);newponto[2] = newpoint(2,0);
			newponto[3] = newnorma(0,0);newponto[4] = newnorma(1,0);newponto[5] = newnorma(2,0);

			MeshAux.push_back(newponto);

			float* ponto = nearpoint(newponto, rootA);

			error += dis(newponto, ponto);

			++j;
		}

		++i;

		error = error/(float) Mesh.size();

		vector<float*> clean1 = vector<float*>();

		Mesh.swap(clean1);

		Mesh = MeshAux;

		cout << "Tentativa: " << i + 1 << endl;
		cout << "Erro: "<< abs(error - preverror) << endl;
		cout << "Distancia media entre as malhas: " << sqrt(dis(centA, centB))<< endl;
		cout << endl;

		if (abs(preverror - error) < threshold){
			target = true;
		}
	}
	delete(rootA);
}
