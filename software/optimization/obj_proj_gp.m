%% Nonlinear variance optimization

function L = obj_proj_gp(q_input, dq_input, uproj, Uref_task, DeltaT, Pj, u_prev, Kv, Gps)
    
    % define variables
    
    var = zeros(1,7);
    
    % projected ddq0
    d2qproj = Pj * (uproj - (Kv * dq_input'));

    % complete ddq
    d2q = Uref_task + d2qproj';

    % Euler integration
    dq = dq_input + DeltaT * d2q;
    q = q_input + DeltaT * dq;   

    % Loss function (to be minimized)
    L = 0.001 * (uproj - u_prev)' * (uproj - u_prev);
    
    for i=1:7
        [~,var(i)] = Gps{i}.predict([q,dq,zeros(1,7)]);
    end

    % update Loss function 
    L = L + 1 * norm(var,2);                % squared norm of variances
%     L = L + max(var);                       % minimize the maximum variance
%     L = L + [0.1 0.1 0.1 0.1 1 1 1]*var';   % Loss with weighted variances
    
end