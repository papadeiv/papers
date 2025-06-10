function f=nodalf(p,u) % SH35 non-linearity F(u) for the 2nd-order system formulation
par=u(p.nu+1:end); 
lam=par(1); 
nup=2;%par(2); 
n=p.nu/2; 
% Define the two new variables (u1,u2) --> (u,Δu)
u1=u(1:n); 
u2=u(n+1:2*n); 
% Define the cubic-quintic non-linearity
f1=(lam-1)*u1 + nup*u1.^3 - u1.^5 - 2*u2; 
f2=0*u2;
% Return the non-lineary F(u)
f=[f1;f2]; 