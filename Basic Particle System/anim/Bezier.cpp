#include <iostream>
#include <vector>
#include <cmath>
#include "BaseSystem.h"
#include <shared/defs.h>
#include <util/util.h>
#include "animTcl.h"
#include <GLmodel/GLmodel.h>
#include "shared/opengl.h"
#include "SampleParticle.h"

using namespace std;


void bezier(double* P0, double* P1, double* P2, double* P3, double t, double *Pt) {
    //This function takes 4 points, a double t and an empty pointer Pt and finds an instantaneous location of Pt given t and the points P0-P3.
	
	Vector Pt1 = (P0 * pow((1 - t), 3));
	Vector Pt2 = (3 * t * pow((1 - t), 2) * P1);
	Vector Pt3 = (3 * pow(t, 2) * (1 - t) * P2);
	Vector Pt4 = (pow(t, 3) * P3);
}
