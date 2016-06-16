function [ R,T ] = ICP( p , q)
    p_mean=mean(p);
    q_mean=mean(q);
    
    p_hat=p-repmat(p_mean,size(p,1),1);
    q_hat=q-repmat(q_mean,size(q,1),1);
    %q_hat=q-q_mean;
    
    M=p_hat'*q_hat;
    [U,~,V] = svd(M);    
    R=U*V';
    T=p_mean'-R*q_mean';
end

