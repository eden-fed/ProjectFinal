function [ gz_enc, gz_bar_enc ] = EvalGzAndGzBar( cageVerteciesB4Map,cageVerteciesAfterMap,NumOfVerticesInEdges )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
deltaS=cageVerteciesB4Map([2:end 1])-cageVerteciesB4Map;
deltaD=cageVerteciesAfterMap([2:end 1])-cageVerteciesAfterMap;

gz = 0.5*(abs(deltaD) + abs(deltaS)) .* deltaD ./ (abs(deltaD).*deltaS); %affine transformation with unit normal
gz_gag = 0.5*(abs(deltaD) - abs(deltaS)) .* deltaD ./ (abs(deltaD).*conj(deltaS)); %affine transformation with unit normal

gz_enc=repelem_ours(gz,NumOfVerticesInEdges);
gz_bar_enc=repelem_ours(gz_gag,NumOfVerticesInEdges);

end

