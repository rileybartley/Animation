#include <GLModel/GLModel.h>
#include <shared/defs.h>
#include <util/util.h>
#include "animTcl.h"
#include "BaseSimulator.h"
#include "ObjectPath.h"

#include <string>
#include <time.h>

// a sample simulator

class PathSimulator : public BaseSimulator
{
public:

	PathSimulator(const std::string& name, ObjectPath* target);
	~PathSimulator();
	
	double calculateRotation( double t );
	double ease(double time, double t);
	int step(double time);
	int init(double time)
	{
		m_object->getState(m_pos0);

		return 0;
	};

	int command(int argc, myCONST_SPEC char** argv) { return TCL_OK; }

protected:

	Vector m_pos0; // initial position
	Vector m_pos;
	double max_speed;
	double current_speed;
	clock_t begin_time = clock();

	

	ObjectPath* m_object;

};