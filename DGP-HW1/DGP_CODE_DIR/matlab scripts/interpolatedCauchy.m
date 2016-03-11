% interpolationCoordinates=inv(interpolationCoordinates);

cvx_begin
    variable f(n,1) complex
    minimize 1
    subjecet to
        A * f == q
cvx_end