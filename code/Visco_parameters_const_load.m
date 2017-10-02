% ----------------------------------------------------------------------
%
%  Dinecsionless formulation and scaling, May 2014, MN
%  Scaled quantities:
% -- charactersitic length L (largest space length)
% -- characteristic coefficient factor M_char = max(E)
% -- characteristic stress             S_shar = 1
% -- characteristic displacement       U_char = L^2/M
% -- characteristic gravity constant   g_char = grav
% -- characteristic density          rho_char = M_char/(g_char*L_char)
% -- characteristic time               T_char = mu_E/eta
%
% Note: rho_char is such that (rho_char*g_char*L_char)/M_char = 1
%       U_char is such that L_char^2/(M_char*U_char) = 1
%       T_char is such that mu_E*T_char/eta = 1
%
% ----------------------------------------------------------------------
% Parameter setting for the visco-elastic problem
% domains   - number of subdomains (no 1 is the one on the surface )
% L         - length (horizontal size) of the domain                (m)
% H         - vertical size of the domain                           (m)
% L_char    - characteristic length of the domain = max(L,H)        (m)
% l_ice     - max width of the ice sheet                            (m)
% h_ice     - max height of the ice sheet                           (m)
% nju       - Poisson ratio (per subdomain)                dimensionless
% E         - Young modulus (per subdomain)      (Pa = N/m^2 = kg/(m.s^2))
% rho_ice   - ice density                                        (kg/m^3)
% rho_earth - Earth density                                      (kg/m^3)
% eta       - viscosity                                          (Pa s)
% grav      - gravity constant                                    (m/s)
% S_char    - characteristic stress (Pa) = 1
% U_char    - characteristic displacement                           (m)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% T_BEG     - time >0, when the ice load starts to be imposed
% T_LGM     - time to reach the glacial maximum
% T_EOG     - time for the ice to melt
% T_run     - time to run after the ice has melted
% load_surf - load surface = 1(boxcar), 2(ellpce), 3(hyperbola)
%
% The array 'Disco' describes the discontinuity regions
% Disco(*,*,k) is the kth region
% Disco(*,1,k) x-interval
% Disco(*,2,k) y-interval
% --------------------------
% Discoef(1,k) = nju    in region k
% Discoef(2,k) = E      in region k
% Discoef(3,k) = mu/eta in region k
%
%    --------------------
%    |                      |    
%    |                      |
%    |   Omega1    |
%    --------------------


%
function [L,H,l_ice,h_ice,rho_ice,rho_earth,...
    Disco,Discoef,grav,load_surf,...
    T_LGM, T_EOG, T, delta_t_char] = Visco_parameters_const_load(domains,wh,Emagn)

global L_char M_char N_char U_char R_char G_char T_char T0
global rhs_const
global advec_const
%% - - - - - - problem and geometry parameters
% - - - - - Tuned to Wu paper "Deformation of an incompressible viscoelastic
% - - - - - - - - flat earth with power law creep" page 37

E_domains = 4.0e+11*ones(domains,1); % Young's modulus, Pa
% E_domains = 1.0e+11*ones(domains,1); % Young's modulus, Pa
E_domains(2,1)=E_domains(2,1)*Emagn;

L0         = 1e+6; %1e+7;   % m
H0(1)      =-5e+5;   % m
H0(2)      =-1e+6;   % m
l_ice0     = 1e+6;   % m
h_ice0     = 2.0e+3;   % m
rho_ice0   =  917;   % kg/m^3
rho_earth0 = 3300;   % kg/m^3
nju(1)     = 0.4999   % dimensionless
nju(2)     = 0.4999;   % dimensionless (to be varied)
eta0       = 1.0e+21;  % Viscosity Pa s was 21
% eta0       = 1.0e+22;  % Viscosity Pa s, Experiment A56
grav0      = 9.81;   % m/s^2
Years      = 18400; % total simulation time (years)
T_LGM0     = 0; %90000*secs_per_year;         % Last Glaciation Maximum duration in seconds.
T_EOG0     = 0; %T_LGM0+10000*secs_per_year;  % End of Glaciation.


load_surf  = []; %18.1e+6;% pa

%% Rescaling
E0            = max(E_domains); % Pa
mju0          = E_domains./(2*(1+nju'));
secs_per_year = 365*24*3600; % Duration of a year in seconds 3.1536e+7
T0            = Years*secs_per_year; % s - total simulation time in seconds

% - - - - - - characteristic values to obtain dimensionless problem
L_char  = L0; %max(abs(L),abs(H));
M_char  = E0;       %max(E_domains) in all subdomains
U_char  = 1;        % not used, dimensionless by def.
G_char  = grav0;
R_char  = rho_earth0;
N_char  = eta0;
% T_char  = 100*secs_per_year; % sinking stops too early
T_char  = 5000*secs_per_year; %

% test for dominating convection:
L_char*grav0*rho_earth0/E0

% - - - - - - scaled values
L         = L0/L_char;
H         = H0/L_char;
E         = E_domains/M_char;
l_ice     = l_ice0/L_char;
h_ice     = h_ice0/L_char;
grav      = grav0/G_char;
rho_ice   = rho_ice0/R_char;
rho_earth = rho_earth0/R_char;
eta       = eta0/N_char; %
T         = ceil(T0/T_char);
T_LGM     = T_LGM0/T_char;
T_EOG     = T_EOG0/T_char;
delta_t_char = secs_per_year/T_char;  % this is one scaled year

Disco(1,1,1:domains) =  0;   %form x
Disco(2,1,1:domains) =  L;   %to   x
%  Horizontal split of the domain into strips (two in this case)
Disco(1,2,1) =  0;     %from y
Disco(2,2,1) = H(1);   %to   y
Disco(1,2,2) = H(1);   %from y
Disco(2,2,2) = H(2);   %to   y

Discoef(1,1:domains) = nju(1:domains);
Discoef(2,1:domains) = E';
mju                  = mju0/M_char;
Maxwell_time_inv     = mju0(1:domains)'/eta0*T_char;
Discoef(3,1:domains) = Maxwell_time_inv;% inverse of the Maxwell time

rhs_const   = L_char/M_char/U_char*G_char*R_char;
advec_const = L_char*G_char*R_char/M_char;
return
