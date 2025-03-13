#include "ObjectPath.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;


ObjectPath::ObjectPath(const std::string& name) :
	BaseSystem(name)
{
	m_sx = 0; 
	m_sy = 0;
	m_sz = 0;
	setVector(m_pos, 0, 0, 0);
	setVector(forward_direction, 0, -1, 0);
	theta = 0.0;

	hermite = new Hermite("hermite");

}	// ObjectPath

void ObjectPath::getState(double* p)
{

	VecCopy(p, m_pos);

}	// ObjectPath::getState

void ObjectPath::getFwdDir(double* p) {
	VecCopy(p, forward_direction);
}

void ObjectPath::setState(double* p)
{

	VecCopy(m_pos, p);

}	// ObjectPath::setState

void ObjectPath::setScale(double sx, double sy, double sz) {
	m_sx = sx;
	m_sy = sy;
	m_sz = sz;
}

void ObjectPath::getPosition(double* p, double t) {
	hermite->getLocation(p, t);
}

int ObjectPath::getLastSampPointIndex() {
	return hermite->lastSampPoint();
}
double ObjectPath::getArcLength(double t) {
	return hermite->arcLength(t);
}

void ObjectPath::reset(double time)
{

	setVector(m_pos, 0, 0, 0);

}	// ObjectPath::Reset


int ObjectPath::command(int argc, myCONST_SPEC char** argv)
{
	if (argc < 1)
	{
		animTcl::OutputMessage("system %s: wrong number of params.", m_name.c_str());
		return TCL_ERROR;
	}
	else if (strcmp(argv[0], "load") == 0)
	{
		//Load points from a file.
		if (argc == 2) {
			
			hermite->loadFromFile(argv[1]);
			Vector origin; 
			hermite -> getPoint(origin, 0);
			setState(origin);
			hermite->populateSampleCurve();
		}
		else
		{
			animTcl::OutputMessage("Usage: load \"<file name>\" ");
			return TCL_ERROR;

		}

	}
	
	glutPostRedisplay();
	return TCL_OK;

}	// ObjectPath::command

void ObjectPath::display(GLenum mode)
{

	glEnable(GL_LIGHTING);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glEnable(GL_COLOR_MATERIAL);

	glTranslated(m_pos[0], m_pos[1], m_pos[2]);
	glRotatef(90, 1, 0, 0);
	glRotatef(theta, 0, 1, 0);
	glScalef(m_sx, m_sy, m_sz);

	if (m_model.numvertices > 0)
		glmDraw(&m_model, GLM_SMOOTH | GLM_MATERIAL);
	else
		glutSolidSphere(1.0, 20, 20);


	glPopMatrix();
	glPopAttrib();

	glPushMatrix();
	glPushAttrib(GL_ALL_ATTRIB_BITS);

	// Draw the points along the curve
	glColor3f(1.0, 0, 0);
	hermite -> displayPoints(4.0);

	// Draw the lines in between the samples
	glColor3f(0.3, 0.7, 0.1);
	hermite ->displayCurve(3.0);
	glPopAttrib();

	
	
	

	

	glPopMatrix();
	glPopAttrib();

}	// ObjectPath::display
