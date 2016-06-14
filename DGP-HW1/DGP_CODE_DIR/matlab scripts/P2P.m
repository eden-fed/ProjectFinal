n=size(C,2); %row length

cvx_begin
    variable  f(n) complex
    minimize(norm((D*f),'fro')) %integral of the norm of the second Derivative
    subject to
        C * f == q
cvx_end