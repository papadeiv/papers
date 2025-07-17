function r=sG(p,u) % RHS G(u) for SH35, see nodalf
f=nodalf(p,u); 
r=p.mat.K*u(1:p.nu)-p.mat.M*f; 