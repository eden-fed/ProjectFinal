function [ R,T ] = CalcTranslationAndRotation4H( p , q)
    p_mean=mean(p);
    q_mean=mean(q);
    
    p_hat=p-repmat(p_mean,size(p,1),1);
    q_hat=q-repmat(q_mean,size(q,1),1);
    
    M=p_hat'*q_hat;
    [U,~,V] = svd(M);    
    R=U*V';
    T=p_mean'-R*q_mean';
    
%     find the angle
%     D=eye(3);
%     D(1:2,1:2)=R;
%     eul=rotm2eul(D);
%     R_angle = eul(1);
end

