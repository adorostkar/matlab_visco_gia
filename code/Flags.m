%Node/edge flags
%
function flag = Flags(type)

Types=strvcat('Dirichlet','Neumann');
flag =strmatch(type,Types);
if flag==0,
   disp('Flags: Wrong type.'),stop
end
return
