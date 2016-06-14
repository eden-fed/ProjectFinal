%if the rank of an NxM matrix is less than min(N,M), then the matrix is singular
R=rank(matToInverse);

if (R >= 3)
	isInvertable=1;
else
	isInvertable=0;
end