% interchange elements of a vector
% E0,E1 are vectors (2,1)
% v - value; fb=1(forward interchange); fb=2 (backward interchange)
%
function E1=my_flip(E0,v,fb)
    
    if E0(fb,1)==v,
        E1=E0;
    else
        E1(1,1)=E0(2,1); E1(2,1)=E0(1,1);
    end
    
    return
