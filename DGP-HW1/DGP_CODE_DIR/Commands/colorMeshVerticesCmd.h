#pragma once

class colorMeshVerticesCmd : public MPxCommand
{
public:
	colorMeshVerticesCmd(); 
	virtual MStatus doIt(const MArgList& argList); 
	static void* creator(); 
	static MSyntax syntax(); //which paramaters are supported
	static MString commandName();

private:
	MStatus ChangeColorByValnce(MFnMesh& meshFn);
	MStatus ChangeColorByCurvature(MFnMesh& meshFn, double min, double max, bool argumentAccepted = false);
	bool checkIfColorSetExists(const MFnMesh& meshFn, const MString& colorSetName);
};
