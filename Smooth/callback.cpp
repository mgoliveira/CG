/*
 * callback.cpp
 *
 *  Created on: 5 de abr de 2016
 *      Author: marco
 */

#include "CHE_L4.hpp"

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

CHE_L4 objeto1;
CHE_L4 objeto2;


void init(){

	glClearColor(0.6, 0.6, 0.6, 1.0);
	glPointSize(5.0);
	glViewport(0, 0, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(dz, 1.5, 0.01, 300);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glEnable(GL_DEPTH_TEST);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	objeto1.read_ply("./PlyFiles/meshclose07.ply"); // Lê a malha a ser suavizada
	objeto2.read_ply("./PlyFiles/meshclose07.ply"); // Armazena a malha original

	//objeto1.write_ply("PlyFiles/cuboantes.ply");

	float Voli = objeto1.compute_volume(); // Calcula o volume inicial da malha

	float n = 1; // número de iterações do laplaciano

	for(int i = 0; i < n; ++i){
		objeto1.laplacianSmooth(Voli, 10.0); // Aplica a suavização. O segundo parâmentro indica o quanto deve suavizar.
		objeto1.compute_COTG();
		objeto1.compute_normals();
	}

	//objeto1.write_ply("PlyFiles/cubodepois.ply");

}

void light(void){

	GLfloat luzAmbiente[4]={0.5,0.5,0.5,1.0};
	GLfloat luzDifusa[4]={1.0,1.0,1.0,1.0};	   // "cor"
	GLfloat luzEspecular[4]={0.6, 0.6, 0.6, 1.0};// "brilho"
	GLfloat posicaoLuz01[4]={3.0, -20.0, -15.0, 1.0};
	GLfloat posicaoLuz02[4]={-3.0, 20.0, -15.0, 1.0};

	// Capacidade de brilho do material
	GLfloat especularidade[4]={0.4,0.4,0.4,1.0};
	GLint especMaterial = 20;

	// Habilita o modelo de colorização de Gouraud
	glShadeModel(GL_SMOOTH);

	// Define a refletância do material
	glMaterialfv(GL_FRONT,GL_SPECULAR, especularidade);
	// Define a concentração do brilho
	glMateriali(GL_FRONT,GL_SHININESS,especMaterial);

	// Ativa o uso da luz ambiente
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, luzAmbiente);

	// Define os parâmetros da luz de número 0
	glLightfv(GL_LIGHT0, GL_AMBIENT, luzAmbiente);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, luzDifusa );
	glLightfv(GL_LIGHT0, GL_SPECULAR, luzEspecular );
	glLightfv(GL_LIGHT0, GL_POSITION, posicaoLuz01 );
	glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.01);


	// Define os parâmetros da luz de número 0
	glLightfv(GL_LIGHT1, GL_AMBIENT, luzAmbiente);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, luzDifusa );
	glLightfv(GL_LIGHT1, GL_SPECULAR, luzEspecular );
	glLightfv(GL_LIGHT1, GL_POSITION, posicaoLuz02 );
	glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, 0.01);

	// Habilita a definição da cor do material a partir da cor corrente
	glEnable(GL_COLOR_MATERIAL);
	//Habilita o uso de iluminação
	glEnable(GL_LIGHTING);
	// Habilita a luz de número 0
	glEnable(GL_LIGHT0);
	// Habilita o depth-buffering
	glEnable(GL_DEPTH_TEST);

	glEnable(GL_POLYGON_OFFSET_FILL);

	glPolygonOffset(1.0, 1.0);

}

void display(void){

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	dx = raio*cos(phi)*sin(theta);
	dy = raio*sin(phi)*sin(theta);
	dz = raio*cos(theta);

	light();

	if(isViewPrincipal == true){
		glLoadIdentity();
		gluLookAt(dx, dy, dz, 0.0, 0.0, 0.0, 0,1,0);
		light();
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

	glColor3f(0.1, 0.4, 0.8);
	glPushMatrix();
	glTranslatef(-10.0, 0.0, 0.0);
	objeto1.draw_smooth();
	glPopMatrix();;
	glColor3f(0.1, 0.4, 0.8);
	glPushMatrix();
	glTranslatef(10.0, 0.0, 0.0);
	objeto2.draw_smooth();
	glPopMatrix();


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
