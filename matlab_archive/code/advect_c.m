% Advection terms element matrix --- buoyancy
% (\div u) c*v
% A(2*nip,2*nip) in separate displacement ordering
% np is the number of points per f.e.
% k  is the current integration point, to evaluate at

function A = advect_c(FUN,FUNT,DER,Coord,np,vec_coeff,nju,wh)
    
    A  =zeros(2*np,2*np);
    if nju==0.5,
        c1 = 0;
        c2 = 0;
    else
        c1 = vec_coeff(3); %c_vec(Coord,1,vec_coeff,wh);
        c2 = vec_coeff(4); %c_vec(Coord,2,vec_coeff,wh);
    end
    
    A(1:np,     1:np)                = c1*FUNT*DER(1,:);
    A(1:np,     np+1:2*np)      = FUNT*(c1*DER(2,:));
    A(np+1:2*np,1:np)           = FUNT*(c2*DER(1,:));
    A(np+1:2*np,np+1:2*np) = c2*FUNT*DER(2,:);
    
    
    return
