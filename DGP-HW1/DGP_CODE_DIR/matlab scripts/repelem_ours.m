function c = repelem_ours(a,b)
a=a';
b=b';
index = zeros(1,sum(b));
index([1 cumsum(b(1:end-1))+1]) = 1;
c = a(cumsum(index));
c=c';
end