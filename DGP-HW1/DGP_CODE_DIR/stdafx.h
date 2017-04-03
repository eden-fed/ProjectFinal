#pragma once


#define short2 MAYA_short2
#define short3 MAYA_short3
#define long2 MAYA_long2
#define long3 MAYA_long3
#define int2 MAYA_int2
#define int3 MAYA_int3
#define float2 MAYA_float2
#define float3 MAYA_float3
#define double2 MAYA_double2
#define double3 MAYA_double3
#define double4 MAYA_double4

//////  MAYA API  //////
#include <maya/MSelectionList.h>
#include <maya/MSyntax.h>
#include <maya/MArgDatabase.h>
#include <maya/MArgList.h>
#include <maya/MPxCommand.h>
#include <maya/MPxHwShaderNode.h>
#include <maya/MMessage.h>
#include <maya/MNodeMessage.h>
#include <maya/MPointArray.h>
#include <maya/MPoint.h>
#include <maya/MBoundingBox.h>
#include <maya/MGlobal.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MImage.h>
#include <maya/MPlugArray.h>
#include <maya/MFnStringData.h>
#include <maya/MEulerRotation.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MMatrix.h>
#include <maya/MPxSurfaceShape.h>
#include <maya/MPxSurfaceShapeUI.h>
#include <maya/MTextureEditorDrawInfo.h>
#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <maya/MItMeshVertex.h>
#include <maya/MItMeshPolygon.h>
#include <maya/MItMeshFaceVertex.h>
#include <maya/MItMeshEdge.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MFloatPointArray.h>
#include <maya/MPxDeformerNode.h>
#include <maya/MItGeometry.h>
#include <maya/MDagModifier.h>
#include <maya/MQuaternion.h>

#undef short2
#undef short3
#undef long2
#undef long3
#undef int2
#undef int3
#undef float2
#undef float3
#undef double2
#undef double3
#undef double4

#define short2 CUDA_short2
#define short3 CUDA_short3
#define long2 CUDA_long2
#define long3 CUDA_long3
#define int2 CUDA_int2
#define int3 CUDA_int3
#define float2 CUDA_float2
#define float3 CUDA_float3
#define double2 CUDA_double2
#define double3 CUDA_double3
#define double4 CUDA_double4

//////  CUDA  //////
#include <cublas.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

#include <helper_cuda.h>
// #include <cutil.h>
// #include <cutil_inline_runtime.h>

#include <cuda_gl_interop.h>

//////  STL  //////
#include <complex>
#include <vector>

//////  Runtime  //////
#include <ctime>
#include <cassert>
#include <limits>


//////  GMM  //////
#include "gmm/gmm.h"

//////  CGAL  //////
 #include <CGAL/Simple_cartesian.h>
 #include <CGAL/Polyhedron_3.h>
 #include <CGAL/Polyhedron_incremental_builder_3.h>
 #include <CGAL/Subdivision_method_3.h>
 #include <CGAL/surface_mesh_parameterization_assertions.h>
 #include <CGAL/Parameterization_mesh_feature_extractor.h>
 #include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/Timer.h>
