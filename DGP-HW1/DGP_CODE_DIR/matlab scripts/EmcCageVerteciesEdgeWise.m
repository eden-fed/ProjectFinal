function [ sampledCageVerteciesAfterMap ] = EmcCageVerteciesEdgeWise( cageVerteciesAfterMap, NumOfVerticesInEdges )

%create q - user sampled cage
%start and end vertecies for each edge
startVertex=cageVerteciesAfterMap;
endVertex=circshift(cageVerteciesAfterMap,size(cageVerteciesAfterMap,1)-1);
%create sampled Cage Vertecies After Map
sampledCageVerteciesAfterMap=zeros(sum(NumOfVerticesInEdges),1);
jj=1;
for ii=1:length(NumOfVerticesInEdges)
    temp=linspace(startVertex(ii),endVertex(ii),NumOfVerticesInEdges(ii)+1);
    temp=temp(1:end-1);
    Len=length(temp);
    sampledCageVerteciesAfterMap(jj:jj+Len-1)=temp;
    jj=jj+Len;
end


end

