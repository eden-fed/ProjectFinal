#pragma once
#include "Utils/GMM_Macros.h"

class reconstructRBFCmd : public MPxCommand
{

public:

	reconstructRBFCmd();
	virtual MStatus	doIt(const MArgList& argList);
	static void* creator();
	static MSyntax syntax();
	static MString commandName();

private:
	unsigned n_flag;
	long grid_flag;
	double eps_flag;

	void setFlags(MArgDatabase& argData, MStatus& stat);
	MStatus getMesh(MArgDatabase& argData, MObject& meshObject);
	GMMDenseColMatrix getDenseMatrixOfDouble(const double value);

};