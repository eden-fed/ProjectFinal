% interpolationCoordinates=inv(interpolationCoordinates);
n = length(q);
cvx_begin
    variable  f(n) complex
    minimize 1
    subject to
        A * f == q
cvx_end