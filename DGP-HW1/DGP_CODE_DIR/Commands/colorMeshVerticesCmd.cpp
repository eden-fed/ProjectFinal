
#include "stdafx.h"

#include "Utils/STL_Macros.h"
#include "Utils/Maya_Macros.h"
#include "Utils/Maya_Utils.h"
//#include "../Utils/MatlabInterface.h"
#include "../Utils/Utilities.h"

#include "colorMeshVerticesCmd.h"

#define VALENCE_COLOR_SET_NAME "Valence"
#define CURVATURE_COLOR_SET_NAME "Curvature"
#define minShortName "-min"
#define minLongName  "-minColor"
#define maxShortName "-max"
#define maxLongName  "-maxColor"
#define PI 3.1416
#define IS_NOT_SET -1

colorMeshVerticesCmd::colorMeshVerticesCmd()
{
}

MStatus colorMeshVerticesCmd::doIt(const MArgList & argList)
{
	MStatus stat = MS::kSuccess;//returned status

	MSyntax commandSyntax = syntax();// check that the syntax is correct
	MArgDatabase argData(commandSyntax, argList, &stat);//if the syntax is correct, insert argList to argData
	MCHECKERROR(stat, "Wrong syntax for command " + commandName());// if it is not correct - error message

	MSelectionList objectsList; //get the list of objects from maya
	stat = argData.getObjects(objectsList); //insert the objects from maya arglist to objectsList
	MCHECKERROR(stat, "Can't access object list");

	MObject object;
	stat = objectsList.getDependNode(0, object);//Get a handle (object) for the 0 indexed element of the selection list.
	MCHECKERROR(stat, "Can't access object");

	MObject meshObject;
	stat = Maya_Utils::getMe_a_Mesh(object, meshObject);//get the mesh model from the object and insert it into meshObject
	MCHECKERROR(stat, "Object is not a mesh");

	MFnMesh meshFn(meshObject, &stat);//get access to meshe's function set
	MCHECKERROR(stat, "Can't access mesh");

	MItMeshPolygon polygonItr(meshFn.object(&stat));//polygon iterator for the mesh
	MCHECKERROR(stat, "failed to instanciate polygon iterator");

	//check that the mesh contains only triangles
	while (!polygonItr.isDone())
	{
		if (3 != polygonItr.polygonVertexCount())//if one of the polygons doesnt have exactly 3 verticies, exit with error
		{
			MCHECKERROR(MS::kFailure, MString("please select mesh only containing triangles. (Found ")
				+ polygonItr.polygonVertexCount() + ")");
			break;
		}
		polygonItr.next();
	}

	//If color sets with the names VALENCE, CURVATURE already exist, their content is eliminated
	meshFn.deleteColorSet(VALENCE_COLOR_SET_NAME); //Deletes a named color set from the object. 

	meshFn.deleteColorSet(CURVATURE_COLOR_SET_NAME);//Deletes a named color set from the object. 

	//create color set with the name VALENCE
	MString Valence = meshFn.createColorSetWithName(VALENCE_COLOR_SET_NAME, (MDGModifier*)0, &stat); //create a color set for valance and get string with its name
	meshFn.setCurrentColorSetName(Valence); //Set the "current" or "working" color set for this object. 
	stat = ChangeColorByValnce(meshFn); 
	MCHECKERROR(stat, "failed to create the Valence color set");

	//create color set with the name CURVATURE
	MString Curvature = meshFn.createColorSetWithName(CURVATURE_COLOR_SET_NAME, (MDGModifier*)0, &stat);//create a color set for curvature and get string with its name
	meshFn.setCurrentColorSetName(Curvature);

	if (argData.isFlagSet(minLongName) && argData.isFlagSet(maxLongName)) //check if the user entered min and max values
	{
		//check the validity of the arguments
		double min, max;
		min = argData.flagArgumentDouble(minLongName, 0, &stat);//Gets the value of the requested flag argument to the given flag as a double.
		MCHECKERROR(stat, "Can't access min arg");
		max = argData.flagArgumentDouble(maxLongName, 0, &stat);//Gets the value of the requested flag argument to the given flag as a double.
		MCHECKERROR(stat, "Can't access max arg");
		if (min > max)
		{
			stat = MS::kFailure;
			MCHECKERROR(stat, "minimum must be below or equal to maximum");
		}
		stat = ChangeColorByCurvature(meshFn, min, max, true);//create the color set by curvature with the inserted values
	}
	else
	stat = ChangeColorByCurvature(meshFn, DBL_MAX, -DBL_MAX); //create the color set by curvature with default values

	MCHECKERROR(stat, "failed to create the curvature color set");

	return stat;
}

void * colorMeshVerticesCmd::creator()
{
	return new colorMeshVerticesCmd;
}


MSyntax colorMeshVerticesCmd::syntax()
{
	MStatus stat = MS::kSuccess;
	MSyntax commandSyntax;

	//flags arguments preceded by a '-' character.
	//command arguments required parameters that follow the flags
	//objects an optional list of Maya objects or the contents of the selection list.

	stat = commandSyntax.addFlag(minShortName, minLongName, MSyntax::kDouble);
	MCHECKERRORNORET(stat, "Can't create Syntax object for this command");

	stat = commandSyntax.addFlag(maxShortName, maxLongName, MSyntax::kDouble);
	MCHECKERRORNORET(stat, "Can't create Syntax object for this command");

	stat = commandSyntax.setObjectType(MSyntax::kSelectionList, 1, 1); //expect exactly one object
	MCHECKERRORNORET(stat, "Can't create Syntax object for this command");

	//useSelectionList 	if true, the command will use the current selection list as arguments if none are specified
	commandSyntax.useSelectionAsDefault(true);

	return commandSyntax;
}

MString colorMeshVerticesCmd::commandName()
{
	return "colorMeshVerticesCmd";
}

//color by the velence of edges
MStatus colorMeshVerticesCmd::ChangeColorByValnce(MFnMesh & meshFn)
{
	MStatus stat = MS::kSuccess;//returned status

	//colorArray - get color by number of edges
	MColor colorCodingDataStructure[7];
	colorCodingDataStructure[0] = MColor(1.f, 0.f, 0.f); 
	colorCodingDataStructure[1] = MColor(0.f, 0.f, 1.f);
	colorCodingDataStructure[2] = MColor(1.f, 1.f, 0.5f);
	colorCodingDataStructure[3] = MColor(0.f, 1.f, 0.f);
	colorCodingDataStructure[4] = MColor(1.f, 0.f, 1.f);
	colorCodingDataStructure[5] = MColor(0.f, 1.f, 1.f);
	colorCodingDataStructure[6] = MColor(0.5f, 0.f, 1.f);

	MColorArray arrVertexColor;	//array with the veticies colors
	arrVertexColor.clear(); 

	int velence = 0; 	// num of neighbors edges
	MIntArray arrVertexIDs;		// array with all the veticies Ids

	MItMeshVertex vertexItr(meshFn.object(&stat));//meshVertex iterator
	MCHECKERROR(stat, "Failed to create iterator");

	while (!vertexItr.isDone()) //loop through all of the verticies
	{
		arrVertexIDs.append(vertexItr.index()); // insert the id of the vertex pointed by the iterator into arrVertexIDs	
		vertexItr.numConnectedEdges(velence); // returns the num of neighbors edges
		if (velence > 9) // if more then 9 -> set to 9
			velence = 9;
		if (velence < 3)
			velence = 3; // if less then 3 -> set to 3
		arrVertexColor.append(colorCodingDataStructure[velence - 3]);// insert the color to the colors array by the number of the neighbors edges
		vertexItr.next();
	}

	stat = meshFn.setVertexColors(arrVertexColor, arrVertexIDs);//match vertex to color
	MCHECKERROR(stat, "faild to set color set ");

	stat = meshFn.setDisplayColors(true);//Determines if the mesh's colors are to be displayed - set to true
	MCHECKERROR(stat, "faild to show colors.");

	stat = meshFn.updateSurface();//Signal that this polygonal mesh has changed and needs to be redrawn.
	MCHECKERROR(stat, "faild to update the surface."); 


	return stat;
}

MStatus colorMeshVerticesCmd::ChangeColorByCurvature(MFnMesh & meshFn, double min, double max, bool argumentAccepted)
{

	MStatus stat = MS::kSuccess;
	MDoubleArray arrCurvaturesOfAllVertecies; //this is an array of curvatures used for color calculation
	MItMeshVertex vertexItr(meshFn.object(&stat)); //meshVertex iterator
	MCHECKERROR(stat, "failed to get iterator for vertices");

	MItMeshPolygon trianglesItr(meshFn.object(&stat)); //MeshPolygon iterator (triangles)
	MCHECKERROR(stat, "failed to get iterator for triangles");

	MIntArray arrVertexIDs; // array with all the veticies Ids
	MIntArray arrNeighborTrianglesToCurrentVertex; // array with all the neighbor triangles Ids
	MIntArray arrTriangleVertices; // array that will hold all the trangle veticies
	MPointArray CrdOfTriangeVertecies;  //this is the coordinates of the trangle veticies

	while (!vertexItr.isDone()) //loop through the trangle veticies and calc the curvature angle
	{
		vertexItr.getConnectedFaces(arrNeighborTrianglesToCurrentVertex); //get neighbor triangles to the vertex point by vertexItr

		double sumOfAngles = 0;
		double vertexCurvature = 0;
		int currentIndex = IS_NOT_SET;

		for (unsigned int j = 0; j < arrNeighborTrianglesToCurrentVertex.length(); ++j)//loop through all the neighbor triangles and add their angle calculation to curve
		{
			int former; //the formar index returned by setIndex function
			stat = trianglesItr.setIndex(arrNeighborTrianglesToCurrentVertex[j], former);// set the index of the current triangle (triangles[j]) to be accessed
			stat = trianglesItr.getVertices(arrTriangleVertices);//get the Id's of the vertices of the current triangle

			trianglesItr.getPoints(CrdOfTriangeVertecies, MSpace::kObject, &stat);//get the positions of the vertices on the current triangle that the iterator is pointing to,  MSpace::kObject is the coordinate system we use.

			//loop through the 3 vertices of the triangle and find the one that coresponds to mesh index vertexItr.index() 
			for (int i = 0; i < 3; ++i)
			{
				if (arrTriangleVertices[i] == vertexItr.index())
				{
					currentIndex = i;
					break;
				}
			}

			//currentIndex was not set than the vertex isnt part of the triangle, issue error
			if (currentIndex == IS_NOT_SET) {
				stat = MS::kFailure;
				MCHECKERROR(stat, "Error finding the vertex inside a face.");
			}

			//set vertex1 and vertex2 to be the other two vertices of the triangle
			int vertex1 = (currentIndex + 1) % 3;
			int vertex2 = (currentIndex + 2) % 3;

			//create two vectors from the vertices to the current vertex
			MVector vector1 = CrdOfTriangeVertecies[vertex1] - CrdOfTriangeVertecies[currentIndex];
			MVector vector2 = CrdOfTriangeVertecies[vertex2] - CrdOfTriangeVertecies[currentIndex];
			vector1.normalize();
			vector2.normalize();
			sumOfAngles += acos(vector1*vector2); //theta=arccos(A*B(Cos(theta))), A and B equal 1
		}

		//check if we are on the boundy of the surface, and compute accordingly
		if (vertexItr.onBoundary()) {
			vertexCurvature = PI - sumOfAngles;
		}
		else {
			vertexCurvature = 2.0f * PI - sumOfAngles;
		}
		
		arrCurvaturesOfAllVertecies.append(sumOfAngles); //add the calculated vertex Curvature to the array of Curvatures
		arrVertexIDs.append(vertexItr.index()); //add the Id of the vertex with that curve to the array of vertices 
		vertexItr.next(); //go to the next vertex
	}

	if (!argumentAccepted) //if we didnt get the min/max arguments, find the min/max values of the angles and coose them.
	{
		for (int i = 0; i < arrCurvaturesOfAllVertecies.length(); i++)
		{
			if (min>arrCurvaturesOfAllVertecies[i])
				min = arrCurvaturesOfAllVertecies[i];
			if (max < arrCurvaturesOfAllVertecies[i]) 
				max = arrCurvaturesOfAllVertecies[i];
		}
	}

	MColorArray colorOfAllVertices;//array with the veticies colors
	colorOfAllVertices.clear();	

	float red, green, blue;
	for (int i = 0; i < arrCurvaturesOfAllVertecies.length(); i++) // get the colors of the curves and insert to colorOfAllVertices
	{
		mapColor(arrCurvaturesOfAllVertecies[i], red, green, blue, min, max);//convert a scalar value (the curve) into R,G,B triplets
		colorOfAllVertices.append(red, green, blue);
	}

	stat = meshFn.setVertexColors(colorOfAllVertices, arrVertexIDs);//Sets the colors of the specified vertices.
	MCHECKERROR(stat, "Failed to configure color set");

	stat = meshFn.setDisplayColors(true);//Determines if the mesh's colors are displayed - set to true
	MCHECKERROR(stat, "Failed to show color set");

	stat = meshFn.updateSurface();//Signal that this polygonal mesh has changed and needs to be redrawn.
	MCHECKERROR(stat, "Failed to updated the surface");

	MString minMax = MString("\n Minimum: ") + min + MString("\n Maximum: ") + max;  //display to Maya’s script editor the min and max values
	MGlobal::displayInfo(minMax);

	return stat;
}