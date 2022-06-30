%% Quadratic problem solved with a nonlinear solver

function L = obj_proj_gp_lin(q_input, dq_input, uproj, Uref_task, DeltaT, Pj, u_prev, Kv, x_lin, Xdata, Ydata, sigma1, sigma2, sigma3, sigma4, sigma5, sigma6, sigma7, beta, sigmaN)
    
    % define variables
    
    sigmaL1 = sigma1(1:21)';
    sigmaL2 = sigma2(1:21)';
    sigmaL3 = sigma3(1:21)';
    sigmaL4 = sigma4(1:21)';
    sigmaL5 = sigma5(1:21)';
    sigmaL6 = sigma6(1:21)';
    sigmaL7 = sigma7(1:21)';
    
    sigmaF1 = sigma1(22);
    sigmaF2 = sigma2(22);
    sigmaF3 = sigma3(22);
    sigmaF4 = sigma4(22);
    sigmaF5 = sigma5(22);
    sigmaF6 = sigma6(22);
    sigmaF7 = sigma7(22);
    
    beta1 = beta(1);
    beta2 = beta(2);
    beta3 = beta(3);
    beta4 = beta(4);
    beta5 = beta(5);
    beta6 = beta(6);
    beta7 = beta(7);
    
    sigmaN1 = sigmaN(1);
    sigmaN2 = sigmaN(2);
    sigmaN3 = sigmaN(3);
    sigmaN4 = sigmaN(4);
    sigmaN5 = sigmaN(5);
    sigmaN6 = sigmaN(6);
    sigmaN7 = sigmaN(7);
    
    var = zeros(1,7);
    
    % projected ddq0
    d2qproj = Pj * (uproj - (Kv * dq_input'));

    % complete ddq
    d2q = Uref_task + d2qproj';

    % Euler integration
    dq = dq_input + DeltaT * d2q;
    q = q_input + DeltaT * dq;   

    % Loss function (to be minimized)
    L = 0.0;
%     L = 0.001 * (uproj - u_prev)' * (uproj - u_prev);

    % Compute linearized variances
    [~, var_lin1] = linearizedRGP(Xdata, x_lin, Ydata(:,1), sigmaL1, sigmaF1, beta1, 1, sigmaN1);
    [~, var_lin2] = linearizedRGP(Xdata, x_lin, Ydata(:,2), sigmaL2, sigmaF2, beta2, 1, sigmaN2);
    [~, var_lin3] = linearizedRGP(Xdata, x_lin, Ydata(:,3), sigmaL3, sigmaF3, beta3, 1, sigmaN3);
    [~, var_lin4] = linearizedRGP(Xdata, x_lin, Ydata(:,4), sigmaL4, sigmaF4, beta4, 1, sigmaN4);
    [~, var_lin5] = linearizedRGP(Xdata, x_lin, Ydata(:,5), sigmaL5, sigmaF5, beta5, 1, sigmaN5);
    [~, var_lin6] = linearizedRGP(Xdata, x_lin, Ydata(:,6), sigmaL6, sigmaF6, beta6, 1, sigmaN6);
    [~, var_lin7] = linearizedRGP(Xdata, x_lin, Ydata(:,7), sigmaL7, sigmaF7, beta7, 1, sigmaN7);

    % compute dx and x_hat
%     dx = [q_input, dq_input, d2q] - x_lin; % delta x (1 by 21 vector)
    
    dx = [q, dq, d2q] - x_lin; % delta x (1 by 21 vector)
    
    x_hat = [1; dx']; % x_hat

    % cut V in order to separate t and t+1
    var_lin1(16:22, 2:15) = zeros(7,14);
    var_lin2(16:22, 2:15) = zeros(7,14);
    var_lin3(16:22, 2:15) = zeros(7,14);
    var_lin4(16:22, 2:15) = zeros(7,14);
    var_lin5(16:22, 2:15) = zeros(7,14);
    var_lin6(16:22, 2:15) = zeros(7,14);
    var_lin7(16:22, 2:15) = zeros(7,14);
    
    
    % compute prediction variances
    var(1) = x_hat' * var_lin1 * x_hat; % first output prediction variance
    var(2) = x_hat' * var_lin2 * x_hat; % second output prediction variance
    var(3) = x_hat' * var_lin3 * x_hat; % third output prediction variance
    var(4) = x_hat' * var_lin4 * x_hat; % first output prediction variance
    var(5) = x_hat' * var_lin5 * x_hat; % second output prediction variance
    var(6) = x_hat' * var_lin6 * x_hat; % third output prediction variance
    var(7) = x_hat' * var_lin7 * x_hat; % first output prediction variance
    
    % update Loss function 
    L = L + 1 * norm(var,2);                % squared norm of variances
%     L = L + max(var);                       % minimize the maximum variance
%     L = L + [0.1 0.1 0.1 0.1 1 1 1]*var';   % Loss with weighted variances
    
end