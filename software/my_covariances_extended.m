%% Compute kernel covariances by matlab api

function [K, K_, K__] = my_covariances_extended(X, x_, sigmaL, sigmaF)
    
    usepdist = 'true';
    dim_input = size(X,2);
    
    % compute euclidean distances
    K = classreg.learning.gputils.calcDistance(X(:,1)/sigmaL(1), X(:,1)/sigmaL(1), usepdist);
    K_ = classreg.learning.gputils.calcDistance(X(:,1)/sigmaL(1), x_(:,1)/sigmaL(1), usepdist);
    K__ = classreg.learning.gputils.calcDistance(x_(:,1)/sigmaL(1), x_(:,1)/sigmaL(1), usepdist);
    for r = 2:dim_input
        K = K + classreg.learning.gputils.calcDistance(X(:,r)/sigmaL(r), X(:,r)/sigmaL(r), usepdist);
        K_ = K_+ classreg.learning.gputils.calcDistance(X(:,r)/sigmaL(r), x_(:,r)/sigmaL(r), usepdist);
        K__ = K__+ classreg.learning.gputils.calcDistance(x_(:,r)/sigmaL(r), x_(:,r)/sigmaL(r), usepdist);    
    end
    
    % apply exponential for computing the kernel
    K = (sigmaF^2)*exp(-0.5*K);
    K_ = (sigmaF^2)*exp(-0.5*K_);
    K__ = (sigmaF^2)*exp(-0.5*K__);
    
end