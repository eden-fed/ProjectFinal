n = length(q);

cvx_begin
    variable  f(n) complex
    minimize 1
    subject to
        C * f == q
cvx_end