#include "CircleSimulator.h"

CircleSimulator::CircleSimulator(const std::string& name, BaseSystem* target) :
	BaseSimulator(name),
	m_object(target)
{
}

CircleSimulator::~CircleSimulator()
{
}


int CircleSimulator::step(double time)
{
	Vector pos;
	m_object->getState(pos);

	double twoPIfreqTime = 2.0 * 3.14 * frequency * time;
	pos[0] = m_pos0[0] + amplitude * cos(twoPIfreqTime);
	pos[1] = m_pos0[1] + amplitude * sin(twoPIfreqTime);

	m_object->setState(pos);

	return 0;

}