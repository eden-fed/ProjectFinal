/////////////////////////////////////
//Copyright (C) Dr. Ofir Weber 2013//
/////////////////////////////////////


#include "stdafx.h"
#include <maya/MFnPlugin.h> //for some reason, I can't put this line in stdafx.h file

#include "Commands/TriangulatePolygonCmd.h"
#include "Commands/colorMeshVerticesCmd.h"
#include "Commands/inverseMatrixCmd.h"
#include "Nodes/SpaceDeformer2D.h"

#include "Utils/Maya_Macros.h"



MStatus initializePlugin(MObject obj)
{ 
	MStatus stat;

	MFnPlugin plugin(obj, "Weber", "2014" , "Any");
	
	REGISTER_NODE(plugin, MPxNode::kDeformerNode, SpaceDeformer2D, NULL);
	REGISTER_COMMAND_WITH_SYNTAX(plugin, TriangulatePolygonCmd);
	REGISTER_COMMAND_WITH_SYNTAX(plugin, colorMeshVerticesCmd);
	REGISTER_COMMAND_WITH_SYNTAX(plugin, inverseMatrixCmd);


	return MS::kSuccess;
}

MStatus uninitializePlugin(MObject obj)
{
	MStatus stat;

	MFnPlugin plugin(obj);

	DEREGISTER_NODE(plugin, SpaceDeformer2D);
	DEREGISTER_COMMAND(plugin, TriangulatePolygonCmd);
	DEREGISTER_COMMAND(plugin, colorMeshVerticesCmd);
	DEREGISTER_COMMAND(plugin, inverseMatrixCmd);


	return MS::kSuccess;
}

