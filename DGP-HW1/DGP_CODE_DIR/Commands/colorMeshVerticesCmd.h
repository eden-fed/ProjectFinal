#pragma once

#define VALENCE_COLOR_SET_NAME "Valence"
#define CURVATURE_COLOR_SET_NAME "Curvature"
#define minShortName "-min"
#define minLongName  "-minColor"
#define maxShortName "-max"
#define maxLongName  "-maxColor"


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
