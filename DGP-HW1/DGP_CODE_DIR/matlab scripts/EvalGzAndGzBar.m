function [ gz, gz_gag ] = EvalGzAndGzBar( cageVerteciesB4Map,cageVerteciesAfterMap,NumOfVerticesInEdges )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
deltaS=cageVerteciesB4Map([2:end 1])-cageVerteciesB4Map;
deltaD=cageVerteciesAfterMap([2:end 1])-cageVerteciesAfterMap;

gz = 0.5*(abs(deltaD) + abs(deltaS)) .* deltaD ./ (abs(deltaD).*deltaS); %affine transformation with unit normal
gz_gag = 0.5*(abs(deltaD) - abs(deltaS)) .* deltaD ./ (abs(deltaD).*conj(deltaS)); %affine transformation with unit normal

if(any(NumOfVerticesInEdges)~=0)
    gz=repelem_ours(gz,NumOfVerticesInEdges);
    gz_gag=repelem_ours(gz_gag,NumOfVerticesInEdges);
end
end

