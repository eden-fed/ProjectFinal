n=size(C,2); %row length

%l=length(increasedVertecies);
%D=ones(l,n);
%iterate over vector increasedVertecies
% i=0;
% for z = increasedVertecies
%     D(i,:)=d(z);
%     i=i+1;
% end


cvx_begin
    variable  f(n) complex
    minimize(norm(D*f),'fro') %integral of the norm of the second Derivative
    subject to
        C * f == q
cvx_end