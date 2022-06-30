function [mu_lin, var_lin] = linearizedRGP(Xdata, x_lin, Ydata, sigmaL, sigmaF, beta, mean_flag, sigma_N)

    % function handles
    my_covariances = @my_covariances_extended;
        
    % compute numeric covariances (on CURRENT x, the one we are linearizing)
    [K, K_, K__] = my_covariances(Xdata, x_lin, sigmaL, sigmaF); 
    
    % adding noise
    K = K + sigma_N^2 * eye(size(K,1));      
    K__ = K__ + sigma_N^2;

    % inverse of K (efficient computation)
    Kinv = K \ eye(size(K,1));  

    % compute all the derivatives
    K_01_Xx = K_.*(Xdata-x_lin)./(sigmaL.^2);
    K_10_xX = - K_' .* (Xdata'-x_lin')./(sigmaL'.^2);
    K_01_xx = zeros(1,21);
    K_10_xx = zeros(21,1);
    K_11_xx = diag(sigmaF^2./sigmaL.^2);
    
    % basis function
    H = ones(size(Xdata,1),1); % with constant basis function
%        H = [ones(size(Xdata,1),1), Xdata]; %with linear basis function
    media = H*beta;

    % compute linearized mean and variance
    mu_lin = [K_';  K_10_xX] * Kinv * (Ydata-(media*mean_flag));
    var_lin = [K__, K_01_xx; K_10_xx, K_11_xx] - [K_'; K_10_xX] * Kinv * [K_, K_01_Xx];
    
end