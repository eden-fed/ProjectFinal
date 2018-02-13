function [  ] = make_mex( )
%  if(exist('csolve', 'file') ~= 3 || exist('localStep', 'file') ~= 3)
%      fullPathName = which('make_csolve.m');
%      [folderName, ~, ~] = fileparts(fullPathName);
%  	cd(folderName);
%  	make_csolve;
%      make_localStep;
%  end
 
 if(exist('treeCumSum', 'file') ~= 3)
     fullPathName = which('treeCumSum.cpp');
     [folderName, ~, ~] = fileparts(fullPathName);
     mex(fullPathName, '-outdir', folderName);
 end
 
end