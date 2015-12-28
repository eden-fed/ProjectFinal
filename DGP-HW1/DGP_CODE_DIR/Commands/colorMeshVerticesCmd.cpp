
#include "stdafx.h"

#include "Utils/STL_Macros.h"
#include "Utils/Maya_Macros.h"
#include "Utils/Maya_Utils.h"
#include "../Utils/MatlabInterface.h"
#include "../Utils/Utilities.h"

#include "colorMeshVerticesCmd.h"

colorMeshVerticesCmd::colorMeshVerticesCmd()
{
}

MStatus colorMeshVerticesCmd::doIt(const MArgList & argList)
{
	MStatus stat = MS::kSuccess;//returned status

	MSyntax commandSyntax = syntax();// chack that the syntax is correct
	MArgDatabase argData(commandSyntax, argList, &stat);//if the syntax is correct, insert argList to argData
	MCHECKERROR(stat, "Wrong syntax for command " + commandName());// if it is not correct - error message

	MSelectionList objectsList; //get the list of objects from maya
	stat = argData.getObjects(objectsList); //insert the objects to objectsList
	MCHECKERROR(stat, "Can't access object list");

	MObject object;
	stat = objectsList.getDependNode(0, object);//get the 0 index object from objectsList
	MCHECKERROR(stat, "Can't access object");

	MObject meshObject;
	stat = Maya_Utils::getMe_a_Mesh(object, meshObject);//insert the mesh model of object to meshObject
	MCHECKERROR(stat, "Object is not a mesh");

	MFnMesh meshFn(meshObject, &stat);//get access to meshe's function set
	MCHECKERROR(stat, "Can't access mesh");

	MItMeshPolygon polygonItr(meshFn.object(&stat));//polygon iterator
	MCHECKERROR(stat, "failed to instanciate polygon iterator");

	//check that the mesh contains only triangles
	while (!polygonItr.isDone())
	{
		if (3 != polygonItr.polygonVertexCount())
		{
			MCHECKERROR(MS::kFailure, MString("please select mesh only containing triangles. (Found ")
				+ polygonItr.polygonVertexCount() + ")");
			break;
		}
		polygonItr.next();
	}

	//If color sets with the names VALENCE, CURVATURE already exist, their content is eliminated
	if (checkIfColorSetExists(meshFn, VALENCE_COLOR_SET_NAME))
		meshFn.deleteColorSet(VALENCE_COLOR_SET_NAME);
	if (checkIfColorSetExists(meshFn, CURVATURE_COLOR_SET_NAME))
		meshFn.deleteColorSet(CURVATURE_COLOR_SET_NAME);

	//create color set with the name VALENCE
	MString Valence = meshFn.createColorSetWithName(VALENCE_COLOR_SET_NAME, (MDGModifier*)0, &stat); //create a color set for valance and get string with its name
	meshFn.setCurrentColorSetName(Valence); 
	stat = ChangeColorByValnce(meshFn); //create the color set by valance
	MCHECKERROR(stat, "failed to create the Valence color set");

	//create color set with the name CURVATURE
	MString Curvature = meshFn.createColorSetWithName(CURVATURE_COLOR_SET_NAME, (MDGModifier*)0, &stat);//create a color set for curvature and get string with its name
	meshFn.setCurrentColorSetName(Curvature);


	stat = ChangeColorByCurvature(meshFn, -2 * 3.1416, 2 * 3.1416); //create the color set by curvature
	MCHECKERROR(stat, "failed to create the curvature color set");



	return stat;
}

void * colorMeshVerticesCmd::creator()
{
	return new colorMeshVerticesCmd;
}

//color by the velence of edges
MStatus colorMeshVerticesCmd::ChangeColorByValnce(MFnMesh & meshFn)
{
	MStatus stat = MS::kSuccess;//returned status

	//colorArray - get color by number of edges
	MColor colorArray[7];
	colorArray[0] = MColor(1.f, 0.f, 0.f); //check 1.0, 0.0, 0.0
	colorArray[1] = MColor(0.f, 0.f, 1.f);
	colorArray[2] = MColor(1.f, 1.f, 0.5f);
	colorArray[3] = MColor(0.f, 1.f, 0.f);
	colorArray[4] = MColor(1.f, 0.f, 1.f);
	colorArray[5] = MColor(0.f, 1.f, 1.f);
	colorArray[6] = MColor(0.5f, 0.f, 1.f);

	MColorArray colorOfAllVertices;	//array with the veticies colors
	colorOfAllVertices.clear(); 
	int EdgesSum = 0; 	// num of neighbors edges
	MIntArray vertices;		// array with all the veticies Ids

	MItMeshVertex vertex(meshFn.object(&stat));//meshVertex iterator
	MCHECKERROR(stat, "Failed to create iterator");

	while (!vertex.isDone())
	{
		vertices.append(vertex.index()); // insert the id of the vertex pointed by the iterator to the new array	
		vertex.numConnectedEdges(EdgesSum); // returns the num of neighbors edges
		if (EdgesSum > 9) // if more then 9 -> set to 9
			EdgesSum = 9;
		if (EdgesSum < 3)
			EdgesSum = 3; // if less then 3 -> set to 3
		colorOfAllVertices.append(colorArray[EdgesSum - 3]);// insert the color to the colors array by the number of the neighbors edges
		vertex.next();
	}

	stat = meshFn.setVertexColors(colorOfAllVertices, vertices);//Sets the colors of the specified vertices.
	MCHECKERROR(stat, "faild to set color set ");

	stat = meshFn.setDisplayColors(true);//Determines if the mesh's colors are displayed - set to true
	MCHECKERROR(stat, "faild to show colors.");

	stat = meshFn.updateSurface();
	MCHECKERROR(stat, "faild to update the surface."); //Signal that this polygonal mesh has changed and needs to be redrawn.


	return stat;
}

//check if a Color Set With The Name colorSetName Exists
bool colorMeshVerticesCmd::checkIfColorSetExists(const MFnMesh& meshFn, const MString& colorSetName)
{
	MStatus stat = MS::kSuccess;
	MStringArray allColorSetNames;

	stat = meshFn.getColorSetNames(allColorSetNames);
	if (stat != MS::kSuccess)
	{
		MGlobal::displayError("Can't retrieve the color set names list");
		return false;
	}

	bool hasThatName = false;
	for (unsigned int i = 0; i < allColorSetNames.length(); i++)
	{
		if (allColorSetNames[i] == colorSetName)
		{
			hasThatName = true;
			break;
		}
	}

	return hasThatName;
}