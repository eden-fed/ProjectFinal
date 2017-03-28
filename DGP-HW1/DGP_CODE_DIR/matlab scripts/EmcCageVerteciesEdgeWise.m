function [ sampledCageVerteciesAfterMap ] = EmcCageVerteciesEdgeWise( cageVerteciesAfterMap, NumOfVerticesInEdges,n )

    if(sum(NumOfVerticesInEdges)==0)%then we have more vertices then the number we want to increase to
        sampledCageVerteciesAfterMap=cageVerteciesAfterMap;
    else
        tic
        %create q - user sampled cage
        %start and end vertecies for each edge
        startVertex=cageVerteciesAfterMap;
        endVertex=cageVerteciesAfterMap([2:end 1]);
        %create sampled Cage Vertecies After Map
        sampledCageVerteciesAfterMap=zeros(n,1);
        jj=1;
        for ii=1:length(NumOfVerticesInEdges)
            temp=linspace(startVertex(ii),endVertex(ii),NumOfVerticesInEdges(ii)+1);
            temp=temp(1:end-1);
            Len=length(temp);
            sampledCageVerteciesAfterMap(jj:jj+Len-1)=temp;
            jj=jj+Len;
        end
    end
end

