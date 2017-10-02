% Computes the current height of the ice cap as a function of time
%

function h = height_ice(wh,x,y,h_ice,time_length,T_BEG,T_LGM,T_EOG)

global test_problem
switch test_problem
    case {0,8}
        h = 0;
    case {9,10}
        h = h_ice; % boxcar with fixed height
% currently all other ice shapes are excluded/not tested
    case 11  % to be developed, ice increase and melt
     T0 = T_LGM-T_BEG;
     T1 = T_EOG-T_LGM;
     if time_length<T_BEG,
       h = 0;
     elseif time_length<T_LGM,
       h  = h_ice/T0*(time_length- T_BEG);
     elseif time_length<T_EOG,
       h  = h_ice/T1*(time_length- T_EOG);
     else
       h = 0;
     end
end

return
