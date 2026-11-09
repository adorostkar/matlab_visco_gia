function [Stress_cur, Stress_E_cur] = SStime_mid_iter(sigOld, sigEOld, dt, U, V, P, Node, Edge, Face, Face_eorder9, Face_Parent, nnode, hx, hy, time, Maxwell_time_inv, k)
    % Compute approximate stress in the middle (center of mass, point 9) of
    % each element. The computation uses another method to get the elastic
    % stress in the currect time and uses the previously computed stresses to
    % construct the new stress.
    
    [~, Stress_E_cur] = SStime_midpoint(U,V,P,Node,Edge,Face,Face_eorder9,Face_Parent,nnode,hx,hy,time,Maxwell_time_inv);
    
    alpha = Maxwell_time_inv;
    
    cm = 1 - alpha*dt*0.5;
    cp = 1 + alpha*dt*0.5;
    c = alpha*dt*0.5;
    
    ea = exp(-alpha*dt);
    if k == 2
        Stress_cur = cm*Stress_E_cur - ea*c*sigEOld;
    else
        Stress_cur = ea*sigOld + cm*Stress_E_cur - ea*cp*sigEOld;
    end
    
    return;