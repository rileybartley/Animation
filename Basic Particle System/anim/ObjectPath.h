#ifndef MY_OBJECTPATH_H
#define MY_OBJECTPATH_H
#define MAX_POINTS 40

#include "BaseSystem.h"
#include "Hermite.h"

#include <string>

#include <assert.h>

/*

base system class. This class is the starting point for every system.

NOTE: This class MUST NEVER be instantiated directly! you must derive your custom systems
	from this one and OVERLOAD the two methods getState and setState

*/

class ObjectPath : public BaseSystem
{
public:

	ObjectPath(const std::string& name);
	virtual void getState(double* p);
	virtual void getFwdDir(double* p);
	virtual void setState(double* p);
	void getPosition(double* p, double t);
	int getLastSampPointIndex();
	double getArcLength(double t);
	void reset(double time);

	void display(GLenum mode = GL_RENDER);

	void readModel(char* fname) { m_model.ReadOBJ(fname); }
	void setTheta(double newTheta) { theta = newTheta;  }
	void setScale(double sx, double sy, double sz);
	void flipNormals(void) { glmReverseWinding(&m_model); }
	int command(int argc, myCONST_SPEC char** argv);

protected:

	//The following variables are used for the model creation and manipulation
	float m_sx;
	float m_sy;
	float m_sz;
	double theta;

	Vector m_pos;
	Vector forward_direction;

	GLMmodel m_model;

	Hermite* hermite;

	//Vector HermiteCurve;
	//int num_points;

};

#endif