% Compute the element mass matrices
% --------------------------------------------------------------------
function M_elem=Assm_ElAdSaddle_quadrTH_Q2Q1_rhs(Gauss_point,Gauss_weight,...
        FUND_all,DERP_all,CoordP)
    
    % global test_problem
    % global lan mju
    
    np    = 4;                    % number of points per f.e.
    npQ2  = 9;                    % number of points per f.e.
    nip   = size(Gauss_point,2);  % nip = number of integration points
    dim   = 2;                    % problem dimension
    ndof  = npQ2*dim;             % ndof = np*dim
    M_elem  = zeros(npQ2,npQ2);   % Mass matrix for the elasticity (one diag block)
    % ------------------------------------------------------------------
    % Gauss_point(2,nip)        Coord(nip,2)
    % FUN(np)                   DER(2,np)          Jac(2,2)     IJac(2,2)
    % Deriv(2,np)               B_cur(np,ndof)     D(np,np)
    %
    for k=1:nip
        FUND   = FUND_all{k};     % (9x1)
        FUNDT  = FUND';
        DERP = DERP_all{k};
        Jac   = DERP*CoordP;
        Det   = det(Jac);
        M_elem9 = FUNDT*FUND;                     % (9,9)=(9,1)*(1,9)
        M_elem = M_elem + Det*Gauss_weight(k)*M_elem9;
    end
    
    return
    
