/*
 * main.cpp
 *
 *  Created on: 5 de abr de 2016
 *      Author: marco
 */



#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>


#include "callback.hpp"

int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowSize(640,480);
    glutInitWindowPosition(100,100);
    glutCreateWindow("Surface Registration");
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutMotionFunc(motion);
    glutMouseFunc(mouse);
    glutReshapeFunc(reshape);

    init();

    glutMainLoop();

    return 0;
}




