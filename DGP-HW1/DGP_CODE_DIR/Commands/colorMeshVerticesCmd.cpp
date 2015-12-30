
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

	if (argData.isFlagSet(minLongName) && argData.isFlagSet(maxLongName)) //check if the user enterded min and max values
	{
		double min, max;
		min = argData.flagArgumentDouble(minLongName, 0, &stat);
		MCHECKERROR(stat, "Can't access min arg");
		max = argData.flagArgumentDouble(maxLongName, 0, &stat);
		MCHECKERROR(stat, "Can't access max arg");
		if (min > max)
		{
			stat = MS::kFailure;
			MCHECKERROR(stat, "minimum must be below or equal to maximum");
		}
		stat = ChangeColorByCurvature(meshFn, min, max, true);//create the color set by curvature with the inserted values
	}
	else
	stat = ChangeColorByCurvature(meshFn, -DBL_MAX, DBL_MAX); //create the color set by curvature with default values

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

	stat = commandSyntax.addFlag(minShortName, minLongName, MSyntax::kDouble);
	MCHECKERRORNORET(stat, "Can't create Syntax object for this command");

	stat = commandSyntax.addFlag(maxShortName, maxLongName, MSyntax::kDouble);
	MCHECKERRORNORET(stat, "Can't create Syntax object for this command");

	stat = commandSyntax.setObjectType(MSyntax::kSelectionList, 1, 1); //expect exactly one object
	MCHECKERRORNORET(stat, "Can't create Syntax object for this command");

	commandSyntax.useSelectionAsDefault(true);

	return commandSyntax;
}

MString colorMeshVerticesCmd::commandName()
{
	return "colorMeshVertices";
}

//color by the velence of edges
MStatus colorMeshVerticesCmd::ChangeColorByValnce(MFnMesh & meshFn)
{
	MStatus stat = MS::kSuccess;//returned status

	//colorArray - get color by number of edges
	MColor colorArray[7];
	colorArray[0] = MColor(1.f, 0.f, 0.f); 
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

	MItMeshVertex vertexItr(meshFn.object(&stat));//meshVertex iterator
	MCHECKERROR(stat, "Failed to create iterator");

	while (!vertexItr.isDone())
	{
		vertices.append(vertexItr.index()); // insert the id of the vertex pointed by the iterator to the new array	
		vertexItr.numConnectedEdges(EdgesSum); // returns the num of neighbors edges
		if (EdgesSum > 9) // if more then 9 -> set to 9
			EdgesSum = 9;
		if (EdgesSum < 3)
			EdgesSum = 3; // if less then 3 -> set to 3
		colorOfAllVertices.append(colorArray[EdgesSum - 3]);// insert the color to the colors array by the number of the neighbors edges
		vertexItr.next();
	}

	stat = meshFn.setVertexColors(colorOfAllVertices, vertices);//Sets the colors of the specified vertices.
	MCHECKERROR(stat, "faild to set color set ");

	stat = meshFn.setDisplayColors(true);//Determines if the mesh's colors are displayed - set to true
	MCHECKERROR(stat, "faild to show colors.");

	stat = meshFn.updateSurface();
	MCHECKERROR(stat, "faild to update the surface."); //Signal that this polygonal mesh has changed and needs to be redrawn.


	return stat;
}

MStatus colorMeshVerticesCmd::ChangeColorByCurvature(MFnMesh & meshFn, double min, double max, bool argumentAccepted)
{

	MStatus stat = MS::kSuccess;
	MDoubleArray curvatures; //this is an array of curvatures used for color calculation
	MItMeshVertex vertexItr(meshFn.object(&stat)); //meshVertex iterator
	MCHECKERROR(stat, "failed to get iterator for vertices");

	MItMeshPolygon trianglesItr(meshFn.object(&stat)); //MeshPolygon iterator (triangles)
	MCHECKERROR(stat, "failed to get iterator for triangles");

	MIntArray vertices; // array with all the veticies Ids
	MIntArray triangles; // array with all the neighbor triangles Ids
	//MIntArray edges; // array with all the neighbor edges Ids
	MIntArray TriangleVertices; // array that will hold all the trangle veticies
	MPointArray TriangeData;  //this is the coordinates of the trangle veticies

	while (!vertexItr.isDone()) //loop through the trangle veticies and calc the curvature angle
	{
		vertexItr.getConnectedFaces(triangles); //get neighbor triangles to the vertex point by vertexItr
	//	vertexItr.getConnectedEdges(edges);		//get neighbor edges to the vertex point by vertexItr

		double curve = 0;
		int currentIndex = -1;

		for (unsigned int j = 0; j < triangles.length(); ++j)//loop through all the neighbor triangles and add their angle calculation to curve
		{
			int former; //the formar index returned by setIndex function
			stat = trianglesItr.setIndex(triangles[j], former);// set the index of the current triangle (triangles[j]) to be accessed
			stat = trianglesItr.getVertices(TriangleVertices);//get the Id's of the vertices of the current triangle

			trianglesItr.getPoints(TriangeData, MSpace::kObject, &stat);//get the positions of the vertices on the current triangle that the iterator is pointing to,  MSpace::kObject is the coordinate system we use.

			//loop through the 3 vertices of the triangle and find the tringel index that coresponds to mash index vertexItr.index() 
			for (int i = 0; i < 3; ++i)
			{
				if (TriangleVertices[i] == vertexItr.index())
				{
					currentIndex = i;
					break;
				}
			}

			if (currentIndex == -1) {
				stat = MS::kFailure;
				MCHECKERROR(stat, "Error finding the vertex inside a face.");
			}

			//set vertex1 and vertex2 to be the other two vertices of the triangle
			int vertex1 = (currentIndex + 1) % 3;
			int vertex2 = (currentIndex + 2) % 3;
			//create two vectors from the vertices to the current vertex
			MVector vector1 = TriangeData[vertex1] - TriangeData[currentIndex];
			MVector vector2 = TriangeData[vertex2] - TriangeData[currentIndex];
			vector1.normalize();
			vector2.normalize();
			curve += acos(vector1*vector2); //get the angle by a dot product of the two vectors and add it to curve
		}

		//now we have the curve of vertexItr
		vertexItr.onBoundary() ? curve = 3.1416 - curve : curve = 2.0f * 3.1416 - curve; //check if we are on the edge of the surface
		
		curvatures.append(curve); //add the calculated curve to the array of curves
		vertices.append(vertexItr.index()); //add the Id of the vertex with that curve to the array of vertices 
		vertexItr.next();
	}

	if (!argumentAccepted) //if we didnt get the min/max arguments, find the min/max values of the angles and coose them.
	{
		for (int i = 0; i < curvatures.length(); i++)
		{
			if (min>curvatures[i]) min = curvatures[i];
			if (max < curvatures[i]) max = curvatures[i];
		}
	}

	MColorArray colorOfAllVertices;//array with the veticies colors
	colorOfAllVertices.clear();	

	float red, green, blue;
	for (int i = 0; i < curvatures.length(); i++) // get the colors of the curves and insert to colorOfAllVertices
	{
		mapColor(curvatures[i], red, green, blue, min, max);//convert a scalar value (the curve) into R,G,B triplets
		colorOfAllVertices.append(red, green, blue);
	}

	stat = meshFn.setVertexColors(colorOfAllVertices, vertices);//Sets the colors of the specified vertices.
	MCHECKERROR(stat, "Failed to configure color set");

	stat = meshFn.setDisplayColors(true);//Determines if the mesh's colors are displayed - set to true
	MCHECKERROR(stat, "Failed to show color set");

	stat = meshFn.updateSurface();//Signal that this polygonal mesh has changed and needs to be redrawn.
	MCHECKERROR(stat, "Failed to updated the surface");

	MString minMax = MString("Minimum:") + min + MString("\n Maximum") + max;  //display to Maya’s script editor the min and max values
	MGlobal::displayInfo(minMax);

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