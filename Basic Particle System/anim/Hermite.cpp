#include "Hermite.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

Hermite::Hermite(const std::string& name) :
	BaseSystem(name)
	
{
	for (int i = 0; i < MAX_POINTS; i++) {
		setVector(points[i], 0, 0, 0);
		setVector(tangents[i], 0, 0, 0);
	}
	for (int j = 0; j < MAX_POINTS * MAX_POINTS; j++) {
		setVector(sampleCurve[j], 0, 0, 0);
		lookupTable[j] = 0.0;
	}
	num_points = 0;
	lookup_index = 0;
	
}	// Hermite

void Hermite::getPoint(double* p, int index)
{

	VecCopy(p, points[index]);

}	// Hermite::getPoint


void Hermite::getTangent(double* p, int index)
{

	VecCopy(p, tangents[index]);

}	// Hermite::getState

int Hermite::lastSampPoint() {
	return lookup_index-1;
}

void Hermite::setPoint(double* p, int index)
{
	if (index >= num_points) {
		addPoint(p);
	}
	points[index][0] = p[0];
	points[index][1] = p[1];
	points[index][2] = p[2];


}	// Hermite::setState


void Hermite::setTangent(double* p, int index)
{
	if (index >= num_points) {
		addTangent(p);
	}
	tangents[index][0] = p[0];
	tangents[index][1] = p[1];
	tangents[index][2] = p[2];


}	// Hermite::setState

void Hermite::cmInit() {
	
	for (int d = 0; d < 3; d++) {
		//Edge cases
		tangents[0][d] = (2.0 * (points[1][d] - points[0][d]) - ((points[2][d] - points[0][d]) / 2.0));
		tangents[num_points - 1][d] = (2.0 * (points[num_points - 2][d] - points[num_points - 3][d]) - ((points[num_points - 1][d] - points[num_points - 3][d]) / 2.0));

		//Itterate from the second point, to the second to last point. 
		for (int i = 1; i < num_points - 1; i++) {
			tangents[i][d] = ((points[i + 1][d] - points[i - 1][d]) / 2.0);
		}
	}
}


double Hermite::evaluateCurveDirect(int d, double t, int i) {
	
	double a = ( -2 * (points[i + 1][d] - points[i][d])) + tangents[i][d] + tangents[i+1][d];

	double b = ( 3 * (points[i + 1][d] - points[i][d])) - (2 * tangents[i][d]) - tangents[i + 1][d];

	double c = tangents[i][d]; 

	double dd = points[i][d];

	return (a * pow(t, 3)) + (b * pow(t, 2)) + (c * t) + dd;
	
}


void Hermite::displayCurve(float r) {
	glPointSize(r);
	Vector curvePoint;
	double stepSize = 1.0 / double(MAX_POINTS);
	glBegin(GL_LINE_STRIP);

	for (int i = 0; i < num_points-1; i++) {
		for (int j = 0; j < MAX_POINTS; j++) {
			double t = j * stepSize;
			
			curvePoint[0] = evaluateCurveDirect(0, t, i);
			curvePoint[1] = evaluateCurveDirect(1, t, i);
			curvePoint[2] = evaluateCurveDirect(2, t, i);

			glVertex3f(curvePoint[0], curvePoint[1], curvePoint[2]);
		}
	}
	glVertex3f(points[num_points - 1][0], points[num_points-1][1], points[num_points - 1][2]);
	glEnd();
}

void Hermite::populateSampleCurve() {
	Vector curvePoint;
	double stepSize = 1.0 / double(MAX_POINTS);

	for (int i = 0; i < num_points - 1; i++) {
		for (int j = 0; j < MAX_POINTS; j++) {
			double t = j * stepSize;

			curvePoint[0] = evaluateCurveDirect(0, t, i);
			curvePoint[1] = evaluateCurveDirect(1, t, i);
			curvePoint[2] = evaluateCurveDirect(2, t, i);

			VecCopy(sampleCurve[num_samples], curvePoint);
			num_samples++;
		}
	}
}

void Hermite::getLocation(double* p, double t) {
	
	//if (num_samples == 0) {
	//	return;
	//}
	//if (t <= 0) return;
	//if (t > 1) return getLocation(p, 1);
	//Retrieve the index value based on the t provided
	int requested_index = (int)(t * ((double)num_samples - 1.0));
	
	VecCopy(p, sampleCurve[requested_index]);

}


void Hermite::displayPoints(float r) {
	glPointSize(r);
	glBegin(GL_POINTS);
	for (int i = 0; i < num_points; i++) {
		glVertex3dv(points[i]);
	}
	glEnd();
}


/*
	NOTE:
	Always call addTangent before addPoint because both functions use num_points 
	but only addPoint() increments num_points. So if you call addPoint first, then you will 
	be inserting 1 index ahead of the proper index for tangents. 

*/

void Hermite::addPoint(double* p)
{
	VecCopy(points[num_points], p);
	num_points++; 

}	// Hermite::addPoint

void Hermite::addTangent(double* p)
{
	VecCopy(tangents[num_points], p);

}	// Hermite::addTangent


void Hermite::reset(double time)
{
	for (int i = 0; i < MAX_POINTS; i++) {
		setVector(points[i], 0, 0, 0);
		setVector(tangents[i], 0, 0, 0);
	}
	num_points = 0; 

}	// Hermite::Reset


double Hermite::arc_distance(double* index1, double* index2) {
	return sqrt(pow(index2[0] - index1[0], 2) + pow(index2[1] - index1[1], 2) + pow(index2[2] - index1[2], 2));
}

double Hermite::arcLength(double t) {
	
	if (lookup_index == 0) {
		populate_lookupTable();
	}
	if (t <= 0) return 0;
	if (t > 1) return arcLength(1);
	//Retrieve the index value based on the t provided
	double entry_dist = 1.0 /(lookup_index-1.0);
	int requested_index = t * (lookup_index-1);
	double U = entry_dist * requested_index;
	double U1 = entry_dist * (requested_index + 1);
	//i is requested_index
	 
	return lookupTable[requested_index] + ((t-U)/(U1-U)*(lookupTable[requested_index+1]- lookupTable[requested_index]));
}

void Hermite::populate_lookupTable(){
	lookup_index = 0;
	double stepSize = 1.0 / MAX_POINTS;
	Vector curCurvePoint;
	Vector nextCurvePoint;

	for (int i = 0; i < num_points-1; i++) {

		for (int j = 0; j < MAX_POINTS-1; j++) {
			double t1 = j * stepSize;
			double t2 = (j + 1) * stepSize;

			curCurvePoint[0] = evaluateCurveDirect(0, t1, i);
			curCurvePoint[1] = evaluateCurveDirect(1, t1, i);
			curCurvePoint[2] = evaluateCurveDirect(2, t1, i);

			nextCurvePoint[0] = evaluateCurveDirect(0, t2, i);
			nextCurvePoint[1] = evaluateCurveDirect(1, t2, i);
			nextCurvePoint[2] = evaluateCurveDirect(2, t2, i);
			
			if (i == 0 && j == 0) {
				lookupTable[lookup_index] = arc_distance(curCurvePoint, nextCurvePoint);
			}
			else {
				lookupTable[lookup_index] = lookupTable[lookup_index-1] + arc_distance(curCurvePoint, nextCurvePoint);
			}
			
			lookup_index++; 
		}
	}
	//HERE
	lookupTable[lookup_index] = lookupTable[lookup_index - 1] + arc_distance(nextCurvePoint, points[num_points-1]);
	lookup_index++;
}

void Hermite::loadFromFile(char* filename) {
	
	fstream file;

	file.open(filename, ios::in);

	if (file.is_open()) {
		string line;
		//Skip the first line
		getline(file, line);
		while (getline(file, line)) {
			char* lineChar = &line[0];
			animTcl::OutputMessage(lineChar);
			stringstream ss(line);

			double x, y, z, sx, sy, sz;
			ss >> x;
			ss >> y;
			ss >> z;

			ss >> sx;
			ss >> sy;
			ss >> sz;
			Vector tangent = { sx,sy,sz };
			addTangent(tangent);

			Vector point = { x,y,z };
			addPoint(point);

		}
		file.close();
	}

}

int Hermite::command(int argc, myCONST_SPEC char** argv)
{

	if (argc < 1)
	{
		animTcl::OutputMessage("system %s: wrong number of params.", m_name.c_str());
		return TCL_ERROR;
	}
	else if (strcmp(argv[0], "set") == 0)
	{
		//system <name> 
		if (argc == 6)
		{
			if (strcmp(argv[1], "point") == 0) {
				//This is where we set the point specified to x,y,z
				Vector tempPoint;
				tempPoint[0] = (double)atof(argv[3]);
				tempPoint[1] = (double)atof(argv[4]);
				tempPoint[2] = (double)atof(argv[5]);

				int index = (int)atoi(argv[2]);
				setPoint(tempPoint, index);
				populate_lookupTable();

			}
			else if (strcmp(argv[1], "tangent") == 0) {
				//This is where we set the tangent specified to x,y,z
				Vector tempTangent;
				tempTangent[0] = (double)atof(argv[3]);
				tempTangent[1] = (double)atof(argv[4]);
				tempTangent[2] = (double)atof(argv[5]);

				int index = (int)atoi(argv[2]);
				setTangent(tempTangent, index);
				populate_lookupTable();
			}
			//set a point/tangent at a provided index to a specified value. 
		}
		else
		{
			animTcl::OutputMessage("Usage: set <point/tangent> <index> <x y z>");
			return TCL_ERROR;
		}
	}
	else if (strcmp(argv[0], "add") == 0)
	{
		if (argc == 8)
		{
			//add point and tangent <x y z sx sy sz> to points and tangents

			Vector tangent = { (double)atof(argv[5]), (double)atof(argv[6]), (double)atof(argv[7]) };
			addTangent(tangent);

			Vector point = { (double)atof(argv[2]), (double)atof(argv[3]), (double)atof(argv[4]) };

			addPoint(point);

		}
		else
		{
			animTcl::OutputMessage("Usage: add point <x y z sx sy sz> ");
			return TCL_ERROR;

		}
	}
	else if (strcmp(argv[0], "cr") == 0)
	{
		if (argc == 1)
		{
			cmInit();
		}
		else
		{
			animTcl::OutputMessage("Usage: cr ");
			return TCL_ERROR;
		}
	}

	else if (strcmp(argv[0], "getArcLength") == 0)
	{
		if (argc == 2)
		{
			// getArcLength <t>
			double t = atof(argv[1]);
			char message[2048];
				
			strcpy(message,"arclength at t: ");
			strcat(message, argv[1]);
			strcat(message, " = ");
			string arcstring = to_string(arcLength(t));
			strcat(message, arcstring.c_str());
			animTcl::OutputMessage(message);
		}
		else
		{
			animTcl::OutputMessage("Usage: getArcLength <t> ");
			return TCL_ERROR;

		}
	}
	else if (strcmp(argv[0], "load") == 0)
	{
		//Load points from a file.
		if (argc == 2) {
			// load "<file name>"
			loadFromFile(argv[1]);

		}
		else
		{
			animTcl::OutputMessage("Usage: load \"<file name>\" ");
			return TCL_ERROR;

		}


	}
	else if (strcmp(argv[0], "export") == 0)
	{
		//Export points to a file.
		if (argc == 2) {
			// export "<file name>"
			ofstream file1;
			file1.open(argv[1]);

			if (file1.is_open())
			{
				file1 << "hermite " << num_points << '\n';
				for (int i = 0; i < num_points; i++) {
					file1 << points[i][0] << ' ' << points[i][1] << ' ' << points[i][2] << ' ';
					file1 << tangents[i][0] << ' ' << tangents[i][1] << ' ' << tangents[i][2];

					if (i != (num_points - 1)) {
						file1 << '\n';
					}

				}

				file1.close();

			}
			else
			{
				animTcl::OutputMessage("Usage: export \"<file name>\" ");
				return TCL_ERROR;

			}
			//Export points to a file
		}

		glutPostRedisplay();
		return TCL_OK;

	}	// Hermite::command
}

void Hermite::display(GLenum mode)
{
	glEnable(GL_LIGHTING);
	glMatrixMode(GL_MODELVIEW);
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glEnable(GL_COLOR_MATERIAL);
	
	// Draw the points along the curve
	//glColor3f(0, 0, 0);
	//displayPoints(3.0);

	// Draw the lines in between the samples
	glColor3f(0.3, 0.7, 0.1);
	displayCurve(3.0);

	glPopAttrib();

}	// Hermite::display
