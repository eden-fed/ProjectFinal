function [ log_fz_i ] = logarithmExtraction( z_i, fz_i, f_i_before, NumOfVerticesInEdges )

z_i_p = z_i([2:end 1]); %z_i_plus_1
e_i = z_i_p - z_i;

e_i_m = e_i([end 1:end-1]);
fz_i_m = fz_i([end 1:end-1]);

f_i=EmcCageVerteciesEdgeWise( f_i_before, NumOfVerticesInEdges );
f_i_p = f_i([2:end 1]);
ef_i = f_i_p - f_i;

%vertex angles in the range [0,2pi]
cornerSourceAngles = angle(e_i./e_i_m) + pi;
cornerTargetAngles = angle(e_i.*fz_i./(e_i_m.*fz_i_m)) + pi;
cornerAnglesChange =  cornerTargetAngles - cornerSourceAngles;

%sanity check

if(abs(sum(cornerAnglesChange)) > 1e-6)
    warning('The turning number of target polygon is invalid!');
end

cornerAnglesChange(1) = 0;
arg_fz_i = angle(ef_i(1)/e_i(1)) + cumsum(cornerAnglesChange); 
ln_abs_fz_i = log(abs(fz_i));
log_fz_i = ln_abs_fz_i + 1i*arg_fz_i;

%sanity check
if(norm(exp(log_fz_i)-exp(log(fz_i)), Inf) > 1e-6)
    warning('The sanity check for log of derivative failed!');
end

end

