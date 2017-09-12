% Advection terms element matrix - pre-stress
% grad(u \cdot b) v 
% A(2*nip,2*nip) in separate displacement ordering
% np is the number of points per f.e.
% k  is the current integration point, to evaluate at

function A = advect_b(FUN,FUNT,DER,Coord,np,vec_coeff,nju,wh)

A  =zeros(2*np,2*np);
b1 = vec_coeff(1); %b_vec(Coord,1,vec_coeff,wh);
b2 = vec_coeff(2); %b_vec(Coord,2,vec_coeff,wh);

A(1:np,     1:np)                = b1*FUNT*DER(1,:);
A(1:np,     np+1:2*np)      = FUNT*(b2*DER(1,:));
A(np+1:2*np,1:np)           = FUNT*(b1*DER(2,:));
A(np+1:2*np,np+1:2*np) = b2*FUNT*DER(2,:);

return
