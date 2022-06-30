%% LQ variance optimization by 'cutting' matrix V for decoupling time t from time t+1

function [H, f] = quadratic_problem(q_input, dq_input, DeltaT, x_lin, Xdata, Ydata, sigma1, sigma2, sigma3, sigma4, sigma5, sigma6, sigma7, beta, sigmaN)

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
    
    N = length(q_input);
    
   % Compute linearized variances
    [~, var_lin1] = linearizedRGP(Xdata, x_lin, Ydata(:,1), sigmaL1, sigmaF1, beta1, 1, sigmaN1);
    [~, var_lin2] = linearizedRGP(Xdata, x_lin, Ydata(:,2), sigmaL2, sigmaF2, beta2, 1, sigmaN2);
    [~, var_lin3] = linearizedRGP(Xdata, x_lin, Ydata(:,3), sigmaL3, sigmaF3, beta3, 1, sigmaN3);
    [~, var_lin4] = linearizedRGP(Xdata, x_lin, Ydata(:,4), sigmaL4, sigmaF4, beta4, 1, sigmaN4);
    [~, var_lin5] = linearizedRGP(Xdata, x_lin, Ydata(:,5), sigmaL5, sigmaF5, beta5, 1, sigmaN5);
    [~, var_lin6] = linearizedRGP(Xdata, x_lin, Ydata(:,6), sigmaL6, sigmaF6, beta6, 1, sigmaN6);
    [~, var_lin7] = linearizedRGP(Xdata, x_lin, Ydata(:,7), sigmaL7, sigmaF7, beta7, 1, sigmaN7);
 
    % formulate the quadratic problem
    b1 = q_input - x_lin(1:7) + dq_input*DeltaT;
    b2 = dq_input - x_lin(8:14);
    b3 = -x_lin(15:21);
    b = [1; b1'; b2'; b3'];
    
    V = var_lin1 + var_lin2 + var_lin3 + var_lin4 + var_lin5 + var_lin6 + var_lin7;

    % cut V in order to avoid t and t+1 variables to be coupled
    V_ = V;
    V_(16:22, 2:15) = zeros(7,14);
    V_(2:15, 16:22) = zeros(14,7);
    
    % weight t and t+1 submatrices
    V_(16:22, 16:22) = 4*V_(16:22, 16:22);
    V_(2:15, 2:15) = V_(2:15, 2:15);


    % w = A*ddq
    A = [zeros(1,N); eye(N)*DeltaT^2 ; eye(N)*DeltaT; eye(N)];
     
    % quadratic problem
    H = A' * (2*V_) * A;
    f = (2*b'*V_*A)';
        
end