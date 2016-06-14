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

% edgeVectors = vertices(endIndices) - vertices(startIndices);
edgeVectors_gpu = (gpuArray(vertices(endIndices)) - gpuArray(vertices(startIndices)))./2;
% integral_on_edges=gather(edgeVectors_gpu);