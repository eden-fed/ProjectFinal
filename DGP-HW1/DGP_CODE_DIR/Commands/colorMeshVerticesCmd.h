#pragma once



class colorMeshVerticesCmd : public MPxCommand
{
public:
	colorMeshVerticesCmd(); 
	virtual MStatus doIt(const MArgList& argList); //This method does the actual work of the command
	static void* creator(); //returns anew instance of the command class
	static MSyntax syntax(); //which paramaters are supported
	static MString commandName(); //returns the command name

private:
	MStatus ChangeColorByValnce(MFnMesh& meshFn);
	MStatus ChangeColorByCurvature(MFnMesh& meshFn, double min, double max);
};
