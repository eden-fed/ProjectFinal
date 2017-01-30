#include "stdafx.h"
#include "Utils/MatlabInterface.h"
#include "Utils/Maya_Macros.h"
#include "Utils/STL_Macros.h"
#include "Utils/Maya_Utils.h"
#include "Utils/Utilities.h"

#include "Utils/MatlabGMMDataExchange.h"
#include "reconstructRBFCmd.h"
#define nShortName "-n"
#define nLongName  "-numOfRBF"
#define epsShortName "-e"
#define epsLongName  "-epsilon"
#define gridShortName "-g"
#define gridLongName  "-grid"

reconstructRBFCmd::reconstructRBFCmd()
{

}

MStatus reconstructRBFCmd::doIt(const MArgList & argList)
{
	MStatus stat = MS::kSuccess;//returned status

	MSyntax commandSyntax = syntax();// check that the syntax is correct
	MArgDatabase argData(commandSyntax, argList, &stat);//if the syntax is correct, insert argList to argData
	MCHECKERROR(stat, "Wrong syntax for command " + commandName());// if it is not correct - error message

	setFlags(argData,stat);//set all flags from argData

	//get the mesh from the selected object
	MObject meshObject;
	stat=getMesh(argData,meshObject);
	MFnMesh meshFn(meshObject, &stat);//get access to meshe's function set
	MCHECKERROR(stat, "Can't access mesh");

	//create vertex and normals arrays for matlab
	MItMeshVertex  vertexItr(meshFn.object(&stat));//vertex iterator for the mesh
	MCHECKERROR(stat, "failed to instanciate polygon iterator");

	int numOfVerteices = meshFn.numVertices(&stat);
	GMMDenseColMatrix verticesCrdForMatlab(numOfVerteices, 3);
	GMMDenseColMatrix normalsCrdForMatlab(numOfVerteices, 3);

	int i = 0;
	do
	{
		MPoint v = vertexItr.position();//get vertex position
		verticesCrdForMatlab(i, 0) = v.x; verticesCrdForMatlab(i, 1) = v.y; verticesCrdForMatlab(i, 2) = v.z;//enter the vertex position to the matlab array

		MVector n;
		vertexItr.getNormal(n);//get normal to vertex
		if (n.length() != 1)
			n.normalize();
		normalsCrdForMatlab(i, 0) = n.x; normalsCrdForMatlab(i, 1) = n.y; normalsCrdForMatlab(i, 2) = n.z;//enter the normal position to the matlab array

		i++;
		vertexItr.next();
	} while (!vertexItr.isDone());

	MatlabInterface::GetEngine().Eval("clear");
	//send all values to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("verticesCrd", verticesCrdForMatlab);
	MatlabGMMDataExchange::SetEngineDenseMatrix("normalsCrd", normalsCrdForMatlab);
	MatlabGMMDataExchange::SetEngineDenseMatrix("n", getDenseMatrixOfDouble(n_flag));
	MatlabGMMDataExchange::SetEngineDenseMatrix("eps", getDenseMatrixOfDouble(eps_flag));
	MatlabGMMDataExchange::SetEngineDenseMatrix("grid", getDenseMatrixOfDouble(grid_flag));

	MatlabInterface::GetEngine().LoadAndRunScriptToString("//madrid.eng.biu.ac.il/e2012/fedidae1/Desktop/ProjectFinal/DGP-HW1/DGP_CODE_DIR/matlab scripts/reconstructRBF.m");
	//MatlabInterface::GetEngine().LoadAndRunScriptToString("C:/Users/Ben-PC/Documents/MySWprojects/DGP/DGP-Eden/matlab scripts/reconstructRBF.m");


	GMMDenseColMatrix faces, vertices, newEpsilon, newN;
	//Get from matlab the vertices and the faces after the matlab script finished his run
	MatlabGMMDataExchange::GetEngineDenseMatrix("faces", faces);
	MatlabGMMDataExchange::GetEngineDenseMatrix("vertices", vertices);
	//Get the new epsilon and num of rbf in case changed by the matlab script
	MatlabGMMDataExchange::GetEngineDenseMatrix("eps", newEpsilon);
	MatlabGMMDataExchange::GetEngineDenseMatrix("n", newN);
	eps_flag = newEpsilon(0, 0);
	n_flag = newN(0, 0);

	int verticesNum = vertices.size() / 3;
	int facesNum = faces.size() / 3;

	MPointArray vertexArray;//array of all the vertices (for the create function)
	for (int i = 0;i < verticesNum;i++)
		vertexArray.append(MPoint(vertices(i,0),vertices(i,1),vertices(i,2)));
	
	MIntArray polygonCounts(facesNum, 3);//array of vertex counts for each face, all initialized to 3 - each face contains 3 vertices
	
	MIntArray polygonConnects(facesNum*3);//array of vertex connections for each polygon
	for (int i = 0;i < facesNum;i++) 
		for (int j = 0;j < 3;j++)
			polygonConnects[3 * i + j] = faces(i, j) - 1;

	MFnMesh reconstructed(meshObject);
	reconstructed.create(verticesNum, facesNum, vertexArray, polygonCounts, polygonConnects);

	//Update the reconstructed mesh to appear in the maya GUI
	MCHECKERROR(reconstructed.updateSurface(), "faild to update the surface.");

	//Print parameters to user
	std::ostringstream outputStream;
	outputStream <<"m = "<< numOfVerteices << " , n = " << n_flag << " , epsilon = " << eps_flag << " , grid resolution = " << grid_flag << endl;
	displayInfo(outputStream.str().c_str());

	displayInfo("Reconstruction completed.");

	return stat;
}

void * reconstructRBFCmd::creator()
{
	return new reconstructRBFCmd;
}

MSyntax reconstructRBFCmd::syntax()
{
	MStatus stat = MS::kSuccess;
	MSyntax commandSyntax;

	//flags arguments preceded by a '-' character.
	//command arguments required parameters that follow the flags
	//objects an optional list of Maya objects or the contents of the selection list.

	stat = commandSyntax.addFlag(nShortName, nLongName, MSyntax::kUnsigned);
	MCHECKERRORNORET(stat, "Can't create Syntax object for this command 1");

	stat = commandSyntax.addFlag(epsShortName, epsLongName, MSyntax::kDouble);
	MCHECKERRORNORET(stat, "Can't create Syntax object for this command 2");

	stat = commandSyntax.addFlag(gridShortName, gridLongName, MSyntax::kLong);
	MCHECKERRORNORET(stat, "Can't create Syntax object for this command 3");

	stat = commandSyntax.setObjectType(MSyntax::kSelectionList, 1, 1); //expect exactly one object
	MCHECKERRORNORET(stat, "Can't create Syntax object for this command 4");

	//useSelectionList 	if true, the command will use the current selection list as arguments if none are specified
	commandSyntax.useSelectionAsDefault(true);

	return commandSyntax;
}

MString reconstructRBFCmd::commandName()
{
	return "reconstructRBFCmd";
}

void reconstructRBFCmd::setFlags(MArgDatabase & argData,MStatus& stat)
{
	n_flag = 1000;
	grid_flag = 50;
	eps_flag = 0;
	//Case where the flag of RBF functions amount was set by user
	if (argData.isFlagSet(nShortName))
	{
		n_flag = argData.flagArgumentDouble(nShortName, 0, &stat);
	}
	//Case where the flag of epsilon size was set by user
	if (argData.isFlagSet(epsShortName))
	{
		eps_flag = argData.flagArgumentDouble(epsShortName, 0, &stat);
	}
	//Case where the flag of grid resolution was set by user
	if (argData.isFlagSet(gridShortName))
	{
		grid_flag = argData.flagArgumentDouble(gridShortName, 0, &stat);
	}
}

MStatus reconstructRBFCmd::getMesh(MArgDatabase & argData, MObject & meshObject)
{
	MStatus stat = MS::kSuccess;//returned status
	MSelectionList objectsList; //get the list of objects from maya
	stat = argData.getObjects(objectsList); //insert the objects from maya arglist to objectsList
	MCHECKERROR(stat, "Can't access object list");

	MObject object;
	stat = objectsList.getDependNode(0, object);//Get a handle (object) for the 0 indexed element of the selection list.
	MCHECKERROR(stat, "Can't access object");

	stat = Maya_Utils::getMe_a_Mesh(object, meshObject);//get the mesh model from the object and insert it into meshObject
	MCHECKERROR(stat, "Object is not a mesh");

	return stat;
}

GMMDenseColMatrix reconstructRBFCmd::getDenseMatrixOfDouble(const double value)
{
	GMMDenseColMatrix valueMatrix(1, 1);
	valueMatrix(0, 0) = value;
	return valueMatrix;
}
