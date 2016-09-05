% clear all
% clc

%we get from the cpp : compInternalVertices(mX1),adjacencyGraph(mXm),rootVertexIndex(z0)

%create a spanning tree to the graph "adjacencyGraph" , starting with "rootVertexIndex"
%disc is a vector of node indices in the order in which they are discovered. 
%pred is a vector of predecessor node indices (listed in the order of the node indices) of the resulting spanning tree. 
%closed is a vector of node indices in the order in which they are closed
[disc, pred, closed] = graphtraverse(adjacencyGraph, rootVertexIndex, 'Directed', false, 'Method', 'BFS');

endIndices = uint32(disc(2:end));
startIndices = uint32(pred(endIndices));
% endIndices_gpu = gpuArray(uint32(disc(2:end)));
% startIndices_gpu = gpuArray(uint32(pred(endIndices_gpu)));

edgeVectors = (vertices(endIndices) - vertices(startIndices))./2;
% edgeVectors_gpu = (gpuArray(vertices(endIndices_gpu)) - gpuArray(vertices(startIndices_gpu)))./2;
