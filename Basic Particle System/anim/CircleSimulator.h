#include <GLModel/GLModel.h>
#include <shared/defs.h>
#include <util/util.h>
#include "animTcl.h"
#include "BaseSimulator.h"
#include "BaseSystem.h"

#include <string>

// a sample simulator

class CircleSimulator : public BaseSimulator
{
public:

	CircleSimulator(const std::string& name, BaseSystem* target);
	~CircleSimulator();

	int step(double time);
	int init(double time)
	{
		m_object->getState(m_pos0);

		// Set this to set the radius of the circle
		amplitude = 2.0;

		//Set this to change the frequency of the circular movement
		frequency = 1.0;

		return 0;
	};

	int command(int argc, myCONST_SPEC char** argv) { return TCL_OK; }

protected:

	Vector m_pos0; // initial position
	Vector m_vel0; // initial velocity
	Vector m_pos;
	Vector m_vel;

	double amplitude;
	double frequency;

	BaseSystem* m_object;

};
