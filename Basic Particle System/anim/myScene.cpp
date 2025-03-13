////////////////////////////////////////////////////
// // Template code for  CSC 473
////////////////////////////////////////////////////

#ifdef WIN32
#include <windows.h>
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <shared/defs.h>

#include "shared/opengl.h"

#include <string.h>
#include <util/util.h>
#include <GLModel/GLModel.h>
#include "anim.h"
#include "animTcl.h"
#include "myScene.h"
#include "ObjectPath.h"
#include "PathSimulator.h"

#include "Hermite.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;


//#include <util/jama/tnt_stopwatch.h>
//#include <util/jama/jama_lu.h>

// register a sample variable with the shell.
// Available types are:
// - TCL_LINK_INT 
// - TCL_LINK_FLOAT

int g_testVariable = 10;

SETVAR myScriptVariables[] = {
	"testVariable", TCL_LINK_INT, (char *) &g_testVariable,
	"",0,(char *) NULL
};


//---------------------------------------------------------------------------------
//			Hooks that are called at appropriate places within anim.cpp
//---------------------------------------------------------------------------------

// start or end interaction
void myMouse(int button, int state, int x, int y)
{

	// let the global resource manager know about the new state of the mouse 
	// button
	GlobalResourceManager::use()->setMouseButtonInfo( button, state );

	if( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN )
	{
		animTcl::OutputMessage(
			"My mouse received a mouse button press event\n");

	}
	if( button == GLUT_LEFT_BUTTON && state == GLUT_UP )
	{
		animTcl::OutputMessage(
			"My mouse received a mouse button release event\n") ;
	}
}	// myMouse

// interaction (mouse motion)
void myMotion(int x, int y)
{

	GLMouseButtonInfo updatedMouseButtonInfo = 
		GlobalResourceManager::use()->getMouseButtonInfo();

	if( updatedMouseButtonInfo.button == GLUT_LEFT_BUTTON )
	{
		animTcl::OutputMessage(
			"My mouse motion callback received a mousemotion event\n") ;
	}

}	// myMotion


void MakeScene(void)
{
	/*
	* bool success;

	// Create system
	Hermite* hermite = new Hermite("hermite");
	//SampleParticle* porsche = new SampleParticle("porsche");

	// Register systems
	success = GlobalResourceManager::use()->addSystem( hermite, true );
	// make sure it was registered successfully
	assert( success );
	*/
	
	
}	// MakeScene

// OpenGL initialization
void myOpenGLInit(void)
{
	animTcl::OutputMessage("Initialization routine was called.");

}	// myOpenGLInit

void myIdleCB(void)
{
	
	return;

}	// myIdleCB

void myKey(unsigned char key, int x, int y)
{
	 animTcl::OutputMessage("My key callback received a key press event\n");
	return;

}	// myKey

int part1(ClientData clientData, Tcl_Interp* interp, int argc, myCONST_SPEC char** argv)
{
	bool success;
	animTcl::OutputMessage("Part 1 Initiated...");
	// Create system
	Hermite* hermite = new Hermite("hermite");
	//SampleParticle* porsche = new SampleParticle("porsche");

	// Register systems
	success = GlobalResourceManager::use()->addSystem(hermite, true);
	// make sure it was registered successfully
	assert(success);

	return TCL_OK;

}

static int part2(ClientData clientData, Tcl_Interp* interp, int argc, myCONST_SPEC char** argv)
{
	bool success;
	animTcl::OutputMessage("Part 2 Initiated...");
	// Create system
	ObjectPath* objectpath = new ObjectPath("objectpath");

	// Register systems
	success = GlobalResourceManager::use()->addSystem(objectpath, true);
	assert(success);

	objectpath->readModel("data/porsche.obj");
	objectpath->setScale(0.01, 0.01, 0.01);

	PathSimulator* pathSim = new PathSimulator("car", objectpath);
	success = GlobalResourceManager::use()->addSimulator(pathSim);
	assert(success);

	return TCL_OK;

}	// testGlobalCommand

void mySetScriptCommands(Tcl_Interp *interp)
{

	// here you can register additional generic (they do not belong to any object) 
	// commands with the shell

	Tcl_CreateCommand(interp, "part1", part1, (ClientData) NULL,
					  (Tcl_CmdDeleteProc *)	NULL);

	Tcl_CreateCommand(interp, "part2", part2, (ClientData)NULL,
		(Tcl_CmdDeleteProc*)NULL);

}	// mySetScriptCommands
