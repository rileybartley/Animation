#ifndef MY_HERMITE_H
#define MY_HERMITE_H
#define MAX_POINTS 40 

/*

	This is a sample system. It accepts the command "read" followed by the
	path to an OBJ model.

*/


#include "BaseSystem.h"
#include <shared/defs.h>
#include <util/util.h>
#include "animTcl.h"
#include <GLmodel/GLmodel.h>
#include <iostream>
#include <fstream>
#include <string>

#include "shared/opengl.h"

// a sample system
class Hermite : public BaseSystem
{

public:
	Hermite(const std::string& name);
	virtual void getPoint(double* p, int index);
	virtual void getTangent(double* p, int index);
	virtual int lastSampPoint();
	virtual void setPoint(double* p, int index);
	virtual void setTangent(double* p, int index);
	virtual void addPoint(double* p);
	virtual void addTangent(double* p);
	
	void cmInit();
	double evaluateCurveDirect(int d, double t, int i);

	void displayCurve(float r);
	void populateSampleCurve();
	void getLocation(double* p, double t);
	void displayPoints(float r);
	

	void reset(double time);

	void display(GLenum mode = GL_RENDER);
	double arc_distance(double* index1, double* index2);
	double arcLength(double t);
	void populate_lookupTable();
	void loadFromFile(char* filename);
	int command(int argc, myCONST_SPEC char** argv);

protected:

	//Set to 40 maximum points, should change out with a global variable MAX_POINTS
	Vector points[MAX_POINTS];
	Vector tangents[MAX_POINTS];
	Vector sampleCurve[MAX_POINTS * MAX_POINTS];

	double lookupTable[MAX_POINTS * MAX_POINTS];
	
	int lookup_index = 0;

	int num_points = 0;
	int num_samples = 0;

};
#endif
#pragma once
