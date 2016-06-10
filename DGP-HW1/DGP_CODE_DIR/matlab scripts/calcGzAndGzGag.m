function [ gz,gz_gag ] = calcGzAndGzGag( deltaS,deltaD )
%deltaS=s2-s1;
%deltaD=d2-d1;


dz=deltaS/abs(deltaS);
dz_vrt=dz*1i;
df=(deltaD)/abs(deltaS);
df_vrt=((deltaD)/abs(deltaD))*1i;

mat=[real(dz)	real(dz_vrt)	0           0;
    0           0               real(dz)	real(dz_vrt);
    imag(dz)    imag(dz_vrt)    0           0;
    0           0               imag(dz)    imag(dz_vrt)];

vec=[real(df); real(df_vrt); imag(df); imag(df_vrt)];

abcd=mat\vec;

gz=0.5*((abcd(1)+abcd(4))+1i*(abcd(3)-abcd(2)));
gz_gag=0.5*((abcd(1)-abcd(4))+1i*(abcd(3)+abcd(2)));
end

