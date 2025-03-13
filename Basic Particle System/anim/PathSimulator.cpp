#include "PathSimulator.h"
#include "util/Quaternion.h"
#include <string>
#include <math.h>
#include <time.h>

using namespace std;

PathSimulator::PathSimulator(const std::string& name, ObjectPath* target) :
	BaseSimulator(name),
	m_object(target)

{
	setVector(m_pos, 0, 0, 0);
	setVector(m_pos0, 0, 0, 0);
	max_speed = 1.111111111;
	current_speed = 0; 
	
}

PathSimulator::~PathSimulator()
{
}

double PathSimulator::calculateRotation( double t ) {
	Vector A;	
	Vector B0;
	Vector B1;
	m_object->getPosition(B0, t);
	m_object->getPosition(B1, t + 0.01);
	m_object->getFwdDir(A);

	Vector B; 
	VecSubtract(B, B1, B0);
	VecNormalize(B);
	double AdotB = VecDotProd(A, B);
	double theta;
	if (B[0] < 0) {
		theta = (-1 * acos(AdotB) * (180 / PI));
		return theta;
	}
	else {
		theta = acos(AdotB) * (180 / PI);
		return theta;
	}
		

}

double PathSimulator::ease(double time, double total_time) {
	
	double t = time/total_time;
	if (1-t < 0.0001) {
		return 1;
	}

	double t1 = 0.1;
	double t2 = 0.9;


	if (t <=  t1) {
		current_speed = max_speed * (t / t1);
		return max_speed * (pow(t, 2) / (2 * t1));
	}
	else if (t >=  0.9) {
		current_speed = max_speed * (1 - ((t - t2) / (1 - t2)));
		return (max_speed * (t1 / 2)) + 
			(max_speed * (t2 - t1)) + 
			((max_speed - (max_speed * ((t - t2) / (1 - t2)) / 2)) * (t - t2));
	}
	else {
		current_speed = max_speed;
		return (max_speed*(t1/2)) + (max_speed*(t-t1));
	}


}

int PathSimulator::step(double time)
{
	double total_distance = m_object->getArcLength(1);
	double distance = total_distance * ease(time, 18);
	Vector newPos;
	
	if (distance < total_distance) {
		double acceptable_error = 0.0001;
		double Xa = 0;
		double Xb = 1;
		double Xc;

		for (int i = 0; i < 20; i++) {
			Xc = (Xa + Xb) / 2;
			if (distance - m_object->getArcLength(Xc) < 0) {
				//Between Xa and Xc
				Xb = Xc;
			}
			else if (distance - m_object->getArcLength(Xc) > 0) {
				//Between Xc and Xb
				Xa = Xc;
			}
			if (acceptable_error > ((Xb - Xa) / (Xa + Xb))) {
				m_object->getPosition(newPos, Xc);
				
				m_object->setState(newPos);
				break;
			}
		}
		m_object->getPosition(newPos, Xc);
		m_object->setState(newPos);
		double rotation_theta = calculateRotation(Xc);
		m_object->setTheta(rotation_theta);
	}
	else {
		m_object->getPosition(newPos, 1);
		m_object->setState(newPos);
	}
	
	if ((float(clock() - begin_time) / CLOCKS_PER_SEC) > 1) {
		animTcl::OutputMessage("Speed of the car is: %f", current_speed);
		begin_time = clock();
	}
	
	return 0;

}