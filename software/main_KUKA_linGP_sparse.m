%% KUKA LWR 4+ Robot - Simulation with uncertain masses by using the linearized Gaussian processes!

% Notes:
% (1) - for real time computations, consider the parallelization for
%       fitting for the 7 GPs
% (2) - for real time computations, consider to use a C++ code 

% QUADPROG OPTIMIZATION, TIME:  154.748467 seconds
% FMINCON OPTIMIZATION, TIME:   4480.769507 seconds

% clear workspace
clear all;
close all;
clc;

% set format
format long;

% add paths
addpath('./Kinematics/');
addpath('./optimization');

%% Utilities

% choose trajectory 
traj = 'heart';  % circular, sinusoidal, heart

% learning utilities
gp_flag = 1;                % {0,1} - set to 1 to use gp learning
opt_flag = 1;               % {0,1} - set to 1 to include null-space optimization
use_quadprog = 1;           % {0,1} - set to 1 to use quadprog function to optimize the variance, 0 to use fmincon

mod_value = 30;
lin_value = 5;              % set when upgrading the linearization point (iterative) 

basis = 'constant';         % set the basis function the GP have to use
mean_flag = 1;              % (0 or 1) - set if you want to compute mu_hat as "...*Y" or "...*(Y-media)" ;

max_dim = 30;               % set the maximum dataset dimension

enable_check_collision = 0;  % set to 1 to check for self-collisions during the motion

% plot utilities
show_standard = 0;  % {0,1} - set to 1 to show the plots of the actual vs desired ee position
show_error = 1;    % {0,1} - set to 1 to show the plots of EE cartesian error on each axis
show_std = 1; % {0,1} - set to 1 to show the plots of prediction standard deviations 
show_scatter = 1;   % {0,1} - set to 1 to show the scatter plot of the actual vs desired ee position
show_animation = 0; % {0,1} - set to 1 to show the animation after the simulation
save_data = 0;
save_workspace = 0;

%% Import robot by using URDF

% real robot
kuka = importrobot('./models/kuka_lwr.urdf'); 

% nominal robot
% kuka_nominal = importrobot('./models/kuka_lwr_controller.urdf');      % uncertainties: +10%, +10%, -10%, -10%, +10%, -10%, +10%
% kuka_nominal = importrobot('./models/kuka_lwr_controller2.urdf');     % uncertainties: -30%, +20%, -20%, -20%, +20%, -20%, +10% 
kuka_nominal = importrobot('./models/kuka_lwr_controller3.urdf');     % uncertainties: -60%, +40%, -50%, -30%, +20%, -20%, +20% 
% kuka_nominal = importrobot('./models/kuka_lwr_controller_circle.urdf');       % uncertainties: still to be computed

% kuka graphics
addVisual(kuka.Base,"Mesh","./visual/base.STL")
addVisual(kuka.Bodies{1},"Mesh","./visual/link_1.STL")
addVisual(kuka.Bodies{2},"Mesh","./visual/link_2.STL")
addVisual(kuka.Bodies{3},"Mesh","./visual/link_3.STL")
addVisual(kuka.Bodies{4},"Mesh","./visual/link_4.STL")
addVisual(kuka.Bodies{5},"Mesh","./visual/link_5.STL")
addVisual(kuka.Bodies{6},"Mesh","./visual/link_6.STL")
addVisual(kuka.Bodies{7},"Mesh","./visual/link_7.STL")

% setup kuka  
kuka.Gravity = [0,0,-9.81];
kuka_nominal.Gravity = [0,0,-9.81];

kuka.DataFormat = 'row';
kuka_nominal.DataFormat = 'row';

%% Task definition

% time & timing law (row vector)
t0 = 0;         % initial time
Ts = 1e-2;      % simulation & control step
T = 5;          % total simulation time
t = t0:Ts:T;    % time
tau_t = t/T;    % normalized time

% cubic timing law, rest-to-rest motion
s = (-2*tau_t.^3 + 3*tau_t.^2);
ds = ((6*t)/T^2 - (6*t.^2)/T^3);
dds = (6/T^2 - (12*t)/T^3);

% desired Cartesian trajectory: 
if strcmp(traj, 'circular') 
    center = [0.5; 0; 0.1];
    radius = 0.15;
    p_ref = center + radius*[cos(4*pi*s); sin(4*pi*s); zeros(size(s))];
    dp_ref = 4*pi*radius*[-sin(4*pi*s).*ds; cos(4*pi*s).*ds; zeros(size(s))];
    ddp_ff = 4*pi*radius*[-4*pi*cos(4*pi*s).*ds.*ds - sin(4*pi*s).*dds; -4*pi*sin(4*pi*s).*ds.*ds + cos(4*pi*s).*dds; zeros(size(s))];
elseif strcmp(traj, 'sinusoidal')
    p_ref = [0.3+0.2*s ; -0.1 + s/3 + 0.2*cos(4*pi*s); 0.2 + 0.25*s];
    dp_ref = [0.2*ds; ds/3-0.2*(4*pi)*sin(4*pi*s).*ds; 0.25*ds];
    ddp_ff = [0.2*dds; dds/3-0.2*(4*pi)^2*cos(4*pi*s).*ds.*ds - 0.2*(4*pi)*sin(4*pi*s).*dds; 0.25*dds];
elseif strcmp(traj, 'heart')
    center = [0.55; 0.1; 0.3];
    p_ref = center  + [(16*sin(2*pi*s).^3)/16; (13*cos(2*pi*s) - 5*cos(2*2*pi*s) - 2*cos(3*2*pi*s) - cos(4*2*pi*s))/16; zeros(size(s)) ]/4;
    
    dp_ref = [(3*ds*pi.*cos(2*pi*s).*sin(2*pi*s).^2)/2;
                ds.*((5*pi*sin(4*pi*s))/16 - (13*pi*sin(2*pi*s))/32 + (3*pi*sin(6*pi*s))/16 + (pi*sin(8*pi*s))/8);
                zeros(size(s)) ];
    
    ddp_ff = [ (3*dds.*ds*pi.*cos(2*pi*s).*sin(2*pi*s).^2)/2 - ds.*(3*ds*pi^2.*sin(2*pi*s).^3 - 6*ds*pi^2.*cos(2*pi*s).^2.*sin(2*pi*s));
                ((5*pi^2*cos(4*pi*s))/4 - (13*pi^2*cos(2*pi*s))/16 + (9*pi^2*cos(6*pi*s))/8 + pi^2*cos(8*pi*s)).*ds.^2 + dds.*((5*pi*sin(4*pi*s))/16 - (13*pi*sin(2*pi*s))/32 + (3*pi*sin(6*pi*s))/16 + (pi*sin(8*pi*s))/8).*ds;
                zeros(size(s))];
end

% inverse kinematics for setting up the first configuration
ik = inverseKinematics('RigidBodyTree', kuka);
weights = [0, 0, 0, 1, 1, 1]; % only position matters
endEffector = 'kuka_lwr_7_link';
qInitial = zeros(1,7); % Use home configuration as the initial guess
point =  p_ref(:,1);
[qSol, solInfo] = ik(endEffector,trvec2tform(point'),weights,qInitial);

%% initializations

% initial configuration - row vector
q0 = qSol;
q = zeros(length(t), 7);
q(1,:) = q0;

% initial velocity (at rest) - row vector
dq0 = zeros(1,7);
dq = zeros(length(t), 7);
dq(1,:) = dq0;

% initial acceleration - row vector
ddq0 = zeros(1,7);
ddq = zeros(length(t)-1, 7);
ddq(1,:) = ddq0;

% Jacobian
J = J_LWR(q(1),q(2),q(3),q(4),q(5),q(6));

% least square damping
damping = 0.01;

% initial FL torque
TauFL = gravityTorque(kuka_nominal,q0);

% gains (cartesian PD)
Kp = diag([100, 100, 100]);
Kd = 2*sqrt(Kp);
Kv = 10;

% collisions
self_collision = [];
clearCollision(kuka.Bodies{7}); % remove collision cylinder on the EE, problems with self collision check!

% null-space optimization settings (nonlinear)
options_qp = optimoptions('quadprog', 'Display', 'off');
options_opt = optimoptions('fmincon','Algorithm','sqp','Display', 'off');
A = [1;-1];
Aproj = blkdiag(A,A,A,A,A,A,A);
bproj = ones(14,1) * 10;
u0proj = zeros(7,1);
DeltaT = 1e-2;

% standard array structures
joints = [q0 dq0]; % contains both the joint configurations and velocities
singular_values = [];
task_vec = [];
accs_ref = [];
accs_redundant = [];
accs = [];
torques = [];
uopt = zeros(7,1);

% learning array structures
predictions = [];
variances = [];
Xdata = [];
Ydata = [];
Xdata_full = [];
Ydata_full = [];
corrective_torques = [];
needed_torques = [];

%% Simulation 
tic
for i = 1:length(t)-1
    disp(i*Ts)
    
    % Build up dataset
    if mod(i, lin_value) ==  0 || i-lin_value == 1
        if i > lin_value
            % sparse gp, keep only the last 'max_dim' points
%             if size(Xdata,1) > max_dim 
%                 Xdata(1:end-max_dim+lin_value,:) = [];
%                 Ydata(1:end-max_dim+lin_value,:) = [];
%             end

            % compute dataset entries
            aux_X = [joints(end-lin_value:end-1,:),accs_ref(end-lin_value+1:end,1:7)];    
            aux_Y = -accs_ref(end-lin_value+1:end,1:7) + accs(end-lin_value+1:end,1:7);    

            % store entries in the dataset
            Xdata = vertcat(Xdata, aux_X);
            Ydata = vertcat(Ydata, aux_Y);
            Xdata_full = vertcat(Xdata_full, aux_X);
            Ydata_full = vertcat(Ydata_full, aux_Y);
        end

        % Tune GP's hyperparams
        if not(isempty(Xdata)) && (gp_flag == 1)
            if mod(i,mod_value) == 0 || i-lin_value == 1  % tune each mod_value step
               
               if size(Xdata, 1) >= max_dim
                   gp_1 = fitrgp(Xdata,Ydata(:,1), 'FitMethod', 'sd', 'PredictMethod', 'sd', 'ActivesetSize', max_dim, 'ActiveSetMethod', 'sgma', 'KernelFunction','ardsquaredexponential');
                   gp_2 = fitrgp(Xdata,Ydata(:,2), 'FitMethod', 'sd', 'PredictMethod', 'sd', 'ActivesetSize', max_dim, 'ActiveSetMethod', 'sgma', 'KernelFunction','ardsquaredexponential');
                   gp_3 = fitrgp(Xdata,Ydata(:,3), 'FitMethod', 'sd', 'PredictMethod', 'sd', 'ActivesetSize', max_dim, 'ActiveSetMethod', 'sgma', 'KernelFunction','ardsquaredexponential');
                   gp_4 = fitrgp(Xdata,Ydata(:,4), 'FitMethod', 'sd', 'PredictMethod', 'sd', 'ActivesetSize', max_dim, 'ActiveSetMethod', 'sgma', 'KernelFunction','ardsquaredexponential');
                   gp_5 = fitrgp(Xdata,Ydata(:,5), 'FitMethod', 'sd', 'PredictMethod', 'sd', 'ActivesetSize', max_dim, 'ActiveSetMethod', 'sgma', 'KernelFunction','ardsquaredexponential');
                   gp_6 = fitrgp(Xdata,Ydata(:,6), 'FitMethod', 'sd', 'PredictMethod', 'sd', 'ActivesetSize', max_dim, 'ActiveSetMethod', 'sgma', 'KernelFunction','ardsquaredexponential');
                   gp_7 = fitrgp(Xdata,Ydata(:,7), 'FitMethod', 'sd', 'PredictMethod', 'sd', 'ActivesetSize', max_dim, 'ActiveSetMethod', 'sgma', 'KernelFunction','ardsquaredexponential');
               else
                   gp_1 = fitrgp(Xdata,Ydata(:,1),'KernelFunction','ardsquaredexponential');
                   gp_2 = fitrgp(Xdata,Ydata(:,2),'KernelFunction','ardsquaredexponential');
                   gp_3 = fitrgp(Xdata,Ydata(:,3),'KernelFunction','ardsquaredexponential');
                   gp_4 = fitrgp(Xdata,Ydata(:,4),'KernelFunction','ardsquaredexponential');
                   gp_5 = fitrgp(Xdata,Ydata(:,5),'KernelFunction','ardsquaredexponential');
                   gp_6 = fitrgp(Xdata,Ydata(:,6),'KernelFunction','ardsquaredexponential');
                   gp_7 = fitrgp(Xdata,Ydata(:,7),'KernelFunction','ardsquaredexponential');
               end
               % select tuned hyperparameters
               if length(gp_1.KernelInformation.KernelParameters) == 22
                   sigma1 = gp_1.KernelInformation.KernelParameters;        
                   sigma2 = gp_2.KernelInformation.KernelParameters;
                   sigma3 = gp_3.KernelInformation.KernelParameters;
                   sigma4 = gp_4.KernelInformation.KernelParameters;        
                   sigma5 = gp_5.KernelInformation.KernelParameters;
                   sigma6 = gp_6.KernelInformation.KernelParameters;
                   sigma7 = gp_7.KernelInformation.KernelParameters;

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

                   beta1 = gp_1.Beta;
                   beta2 = gp_2.Beta;
                   beta3 = gp_3.Beta;
                   beta4 = gp_4.Beta;
                   beta5 = gp_5.Beta;
                   beta6 = gp_6.Beta;
                   beta7 = gp_7.Beta;

                   sigmaN1 = gp_1.Sigma;
                   sigmaN2 = gp_2.Sigma;
                   sigmaN3 = gp_3.Sigma;
                   sigmaN4 = gp_4.Sigma;
                   sigmaN5 = gp_5.Sigma;
                   sigmaN6 = gp_6.Sigma;
                   sigmaN7 = gp_7.Sigma;
                end 
            else
               if i > max_dim
                   gp_1 = fitrgp(Xdata,Ydata(:,1),'ActiveSetSize', max_dim, 'ActiveSetMethod', 'sgma', 'FitMethod', 'None', 'PredictMethod', 'sd','KernelFunction','ardsquaredexponential','KernelParameters',[sigmaL1';sigmaF1], 'BasisFunction', basis, 'Beta', beta1, 'Sigma', sigmaN1);
                   gp_2 = fitrgp(Xdata,Ydata(:,2),'ActiveSetSize', max_dim, 'ActiveSetMethod', 'sgma', 'FitMethod', 'None', 'PredictMethod', 'sd','KernelFunction','ardsquaredexponential','KernelParameters',[sigmaL2';sigmaF2], 'BasisFunction', basis, 'Beta', beta2, 'Sigma', sigmaN2);
                   gp_3 = fitrgp(Xdata,Ydata(:,3),'ActiveSetSize', max_dim, 'ActiveSetMethod', 'sgma', 'FitMethod', 'None', 'PredictMethod', 'sd','KernelFunction','ardsquaredexponential','KernelParameters',[sigmaL3';sigmaF3], 'BasisFunction', basis, 'Beta', beta3, 'Sigma', sigmaN3);
                   gp_4 = fitrgp(Xdata,Ydata(:,4),'ActiveSetSize', max_dim, 'ActiveSetMethod', 'sgma', 'FitMethod', 'None', 'PredictMethod', 'sd','KernelFunction','ardsquaredexponential','KernelParameters',[sigmaL4';sigmaF4], 'BasisFunction', basis, 'Beta', beta4, 'Sigma', sigmaN4);
                   gp_5 = fitrgp(Xdata,Ydata(:,5),'ActiveSetSize', max_dim, 'ActiveSetMethod', 'sgma', 'FitMethod', 'None', 'PredictMethod', 'sd','KernelFunction','ardsquaredexponential','KernelParameters',[sigmaL5';sigmaF5], 'BasisFunction', basis, 'Beta', beta5, 'Sigma', sigmaN5);
                   gp_6 = fitrgp(Xdata,Ydata(:,6),'ActiveSetSize', max_dim, 'ActiveSetMethod', 'sgma', 'FitMethod', 'None', 'PredictMethod', 'sd','KernelFunction','ardsquaredexponential','KernelParameters',[sigmaL6';sigmaF6], 'BasisFunction', basis, 'Beta', beta6, 'Sigma', sigmaN6);
                   gp_7 = fitrgp(Xdata,Ydata(:,7),'ActiveSetSize', max_dim, 'ActiveSetMethod', 'sgma', 'FitMethod', 'None', 'PredictMethod', 'sd','KernelFunction','ardsquaredexponential','KernelParameters',[sigmaL7';sigmaF7], 'BasisFunction', basis, 'Beta', beta7, 'Sigma', sigmaN7);
               else
                   gp_1 = fitrgp(Xdata,Ydata(:,1),'ActiveSetSize', i, 'ActiveSetMethod', 'sgma', 'FitMethod', 'None', 'PredictMethod', 'sd','KernelFunction','ardsquaredexponential','KernelParameters',[sigmaL1';sigmaF1], 'BasisFunction', basis, 'Beta', beta1, 'Sigma', sigmaN1);
                   gp_2 = fitrgp(Xdata,Ydata(:,2),'ActiveSetSize', i, 'ActiveSetMethod', 'sgma', 'FitMethod', 'None', 'PredictMethod', 'sd','KernelFunction','ardsquaredexponential','KernelParameters',[sigmaL2';sigmaF2], 'BasisFunction', basis, 'Beta', beta2, 'Sigma', sigmaN2);
                   gp_3 = fitrgp(Xdata,Ydata(:,3),'ActiveSetSize', i, 'ActiveSetMethod', 'sgma', 'FitMethod', 'None', 'PredictMethod', 'sd','KernelFunction','ardsquaredexponential','KernelParameters',[sigmaL3';sigmaF3], 'BasisFunction', basis, 'Beta', beta3, 'Sigma', sigmaN3);
                   gp_4 = fitrgp(Xdata,Ydata(:,4),'ActiveSetSize', i, 'ActiveSetMethod', 'sgma', 'FitMethod', 'None', 'PredictMethod', 'sd','KernelFunction','ardsquaredexponential','KernelParameters',[sigmaL4';sigmaF4], 'BasisFunction', basis, 'Beta', beta4, 'Sigma', sigmaN4);
                   gp_5 = fitrgp(Xdata,Ydata(:,5),'ActiveSetSize', i, 'ActiveSetMethod', 'sgma', 'FitMethod', 'None', 'PredictMethod', 'sd','KernelFunction','ardsquaredexponential','KernelParameters',[sigmaL5';sigmaF5], 'BasisFunction', basis, 'Beta', beta5, 'Sigma', sigmaN5);
                   gp_6 = fitrgp(Xdata,Ydata(:,6),'ActiveSetSize', i, 'ActiveSetMethod', 'sgma', 'FitMethod', 'None', 'PredictMethod', 'sd','KernelFunction','ardsquaredexponential','KernelParameters',[sigmaL6';sigmaF6], 'BasisFunction', basis, 'Beta', beta6, 'Sigma', sigmaN6);
                   gp_7 = fitrgp(Xdata,Ydata(:,7),'ActiveSetSize', i, 'ActiveSetMethod', 'sgma', 'FitMethod', 'None', 'PredictMethod', 'sd','KernelFunction','ardsquaredexponential','KernelParameters',[sigmaL7';sigmaF7], 'BasisFunction', basis, 'Beta', beta7, 'Sigma', sigmaN7);
               end
            end            
            Gps = {gp_1,gp_2,gp_3,gp_4,gp_5,gp_6,gp_7}; 
        end 
    end
    
    % compute actual ee position 
    p = f(q(i,1), q(i,2), q(i,3), q(i,4), q(i,5), q(i,6));
    
    % compute jacobian and its derivative
    Jold = J;
    J = J_LWR(q(i,1), q(i,2), q(i,3), q(i,4), q(i,5), q(i,6));
    dJ = dJdt_LWR(dq(i,1), dq(i,2), dq(i,3), dq(i,4), dq(i,5), dq(i,6), q(i,1), q(i,2), q(i,3), q(i,4), q(i,5), q(i,6));
    
    % check for singular values
    [~,S,~] = svd(J);
    singular_value = min(diag(S));
    
    % actual ee velocity and acceleration
    dp = J * dq(i,:)';
    ddp = dJ * dq(i,:)' + J * ddq(i,:)';
    
    % cartesian errors
    err_p = p_ref(:,i) - p;
    err_dp = dp_ref(:,i) - dp;
    
    % acceleration reference
    ddp_ref = Kp * err_p + Kd * err_dp + ddp_ff(:,i);
    
    % projector
    Pj = eye(7) - pinv(J) * J;
    
    % damp task accelerations    
    A = J' * J + damping * eye(length(q(i,:)));
    B = J' * (ddp_ref - dJ * dq(i,:)');              
    X1 = lsqminnorm(A,B);    
    Uref_task = X1';
    
    %primo ordine
    % J dq = dr
    % dq = pinv(J) dr
    % pinv  -> min norm2 J dq - dr 
    % dls   -> min (norm2 J dq - dr + mu norm2 dq)
%           -> soluzione = J
    
    % secondo ordine
    % ddq = pinv(J) ddx = pinv(J) (ddr - dJdq)
    
    % dls -> min (norm 2 J ddq - ddx + mu norm2 ddq)
%             -> soluzione:  
    
    
    
    % (J_damp)y - J'(ddx) -> y = inv(J'J) (J' ddx)
    % J_DLS ddr - ddq = 0
    % 
    % preferred acceleration - null-space optimization 
    
%   ---------------------------------------------------------------------------------------
    if((gp_flag == 1) && (opt_flag == 1) && not(isempty(Xdata)) && (length(gp_1.KernelInformation.KernelParameters) == 22))
        
        % set inputs to give to fmincon/quadprog
        beta = [beta1, beta2, beta3, beta4, beta5, beta6, beta7];
        sigmaN = [sigmaN1, sigmaN2, sigmaN3, sigmaN4, sigmaN5, sigmaN6, sigmaN7];
        
        % set linearization point
        if (mod(i,lin_value) == 0 || i-lin_value == 1)
            x_lin = Xdata(end,:);
        end
        
        if (use_quadprog == 0)  % use fmincon 
            u_prev = uopt;
            [uopt,fval,exitflag,output] = fmincon(@(uproj)obj_proj_gp_lin(q(i,:), dq(i,:), uproj, Uref_task, DeltaT, Pj, u_prev, Kv, x_lin, Xdata, Ydata, sigma1, sigma2, sigma3, sigma4, sigma5, sigma6, sigma7, beta, sigmaN), u0proj, Aproj, bproj, [], [], [], [],[], options_opt);
            u0proj = uopt;
        else                    % use quadprog
            [H, f_] = quadratic_problem_2(q(i,:), dq(i,:), DeltaT, x_lin, Xdata, Ydata, sigma1, sigma2, sigma3, sigma4, sigma5, sigma6, sigma7, beta, sigmaN);
            delta_ddq = quadprog((H+H')/2, f_, Aproj, bproj, [], [], [], [], [], options_qp);
            uopt = delta_ddq;
        end
        
    else 
        uopt = zeros(7,1);
    end
%   ---------------------------------------------------------------------------------------
    
    
    % redundancy resolution: null-space optimization
    Uref_redundant = (Pj * (uopt - Kv*dq(i,:)') )';
    
    % total reference acceleration
    Uref = Uref_task + Uref_redundant;
    
    % GP predictions
    if(not(isempty(Xdata)) && (gp_flag == 1) && length(gp_1.KernelInformation.KernelParameters) == 22 )    
        
        % upgrade and set linearization point
        if (mod(i,lin_value) == 0 || i-lin_value == 1)
            x_lin = Xdata(end,:);
        end

        % call the linearized RGP
        [mu_lin1, var_lin1] = linearizedRGP_sparse(gp_1.ActiveSetVectors, x_lin, gp_1.Alpha, sigmaL1, sigmaF1, beta1, sigmaN1);
        [mu_lin2, var_lin2] = linearizedRGP_sparse(gp_2.ActiveSetVectors, x_lin, gp_2.Alpha, sigmaL2, sigmaF2, beta2, sigmaN2);
        [mu_lin3, var_lin3] = linearizedRGP_sparse(gp_3.ActiveSetVectors, x_lin, gp_3.Alpha, sigmaL3, sigmaF3, beta3, sigmaN3);
        [mu_lin4, var_lin4] = linearizedRGP_sparse(gp_4.ActiveSetVectors, x_lin, gp_4.Alpha, sigmaL4, sigmaF4, beta4, sigmaN4);
        [mu_lin5, var_lin5] = linearizedRGP_sparse(gp_5.ActiveSetVectors, x_lin, gp_5.Alpha, sigmaL5, sigmaF5, beta5, sigmaN5);
        [mu_lin6, var_lin6] = linearizedRGP_sparse(gp_6.ActiveSetVectors, x_lin, gp_6.Alpha, sigmaL6, sigmaF6, beta6, sigmaN6);
        [mu_lin7, var_lin7] = linearizedRGP_sparse(gp_7.ActiveSetVectors, x_lin, gp_7.Alpha, sigmaL7, sigmaF7, beta7, sigmaN7);

        % compute dx and x_hat
        dx = [q(i,:), dq(i,:), Uref] - x_lin; % delta x (1 by 21 vector)
        x_hat = [1; dx']; % x_hat

        % compute predictions mean 
        new_mu1 = mu_lin1'*x_hat + beta1;  
        new_mu2 = mu_lin2'*x_hat + beta2;  
        new_mu3 = mu_lin3'*x_hat + beta3;  
        new_mu4 = mu_lin4'*x_hat + beta4;  
        new_mu5 = mu_lin5'*x_hat + beta5;  
        new_mu6 = mu_lin6'*x_hat + beta6;  
        new_mu7 = mu_lin7'*x_hat + beta7;  

        % compute predictions variance
        new_var1 = x_hat' * var_lin1 * x_hat;
        new_var2 = x_hat' * var_lin2 * x_hat; 
        new_var3 = x_hat' * var_lin3 * x_hat; 
        new_var4 = x_hat' * var_lin4 * x_hat;
        new_var5 = x_hat' * var_lin5 * x_hat; 
        new_var6 = x_hat' * var_lin6 * x_hat; 
        new_var7 = x_hat' * var_lin7 * x_hat;
    
        % regroup predictions and variances
        prediction_vector = [new_mu1, new_mu2, new_mu3, new_mu4, new_mu5, new_mu6, new_mu7];
        variance_vector = [new_var1, new_var2, new_var3, new_var4, new_var5, new_var6, new_var7];
        
        % subtract the GP prediction to the reference acceleration
        Uref = Uref - prediction_vector;
    else
        prediction_vector = zeros(1,7);
        variance_vector  = zeros(1,7);
    end
    
    % check for self collisions
    if enable_check_collision == 1 && i > 1
        self_collision(i) = checkCollision(kuka, q(i,:));
    end
    check = any(self_collision == 1);
    
    if check == 0
        % Nominal control law: Computed Torque (FL + PD + FFW) - row vector
        TauFL = (massMatrix(kuka_nominal, q(i,:)) * Uref')' + velocityProduct(kuka_nominal,q(i,:),dq(i,:)) + gravityTorque(kuka_nominal,q(i,:));
    else
        % Collision detected!
        disp('Self collision detected!')
        break;
    end
    
    % dynamic model inversion - row vector
    ddq(i,:) = ( inv(massMatrix(kuka,q(i,:))) * (TauFL - velocityProduct(kuka,q(i,:),dq(i,:)) - gravityTorque(kuka,q(i,:)))' )';
  
    % numeric integration (Euler method)
    dq(i+1,:) = dq(i,:) + ddq(i,:) * Ts;      
    q(i+1,:) = q(i,:) + dq(i,:) * Ts;
       
    % Compute DeltaM and DeltaNi 
    DeltaM = massMatrix(kuka,q(i,:)) - massMatrix(kuka_nominal,q(i,:));            
    DeltaNi =  velocityProduct(kuka,q(i,:),dq(i,:)) + gravityTorque(kuka,q(i,:)) - velocityProduct(kuka_nominal,q(i,:),dq(i,:)) - gravityTorque(kuka_nominal,q(i,:)); 
    
    % Corrective Torque
    Torque_corrective = (-massMatrix(kuka_nominal, q(i,:)) * prediction_vector')';
    
    % Real torque disturbance acting on the system
    Torque_needed = (DeltaM * ddq(i,:)' + DeltaNi')';
    
    % update standard array structures
    task_vec = vertcat(task_vec, [p',dp',ddp',p_ref(:,i)',dp_ref(:,i)',ddp_ref']);
    accs_ref = vertcat(accs_ref, Uref);
    accs_redundant = vertcat(accs_redundant, uopt');
    accs = vertcat(accs, ddq(i,:));
    singular_values = vertcat(singular_values, singular_value);
    joints = vertcat(joints, [q(i+1,:), dq(i+1,:)]);
    
    % update learning array structures
    predictions = vertcat(predictions, prediction_vector);
    variances = vertcat(variances, variance_vector);
    corrective_torques = vertcat(corrective_torques, Torque_corrective);
    needed_torques = vertcat(needed_torques, Torque_needed);
    torques = vertcat(torques,TauFL);
    
end
toc 

time = toc;

%% plot variables
if check == 1
    
    show_standard = 0;
    show_error = 0;
    show_std = 0;
    show_scatter = 0;
    show_animation = 0;
    
    path = p_ref';
    
    % show the self-colliding configuration
    collisionIndex = find(self_collision, 1);
    show(kuka, q(collisionIndex,:));
    view(3)
    ax = gca;
    ax.Projection = 'orthographic';
    hold on
    plot3(path(:,1), path(:,2), path(:,3), 'r','Linewidth',2)
    title('Self colliding configuration')
end

%% Results (errors and torques)

if show_standard == 1
    
    % EE position x-axis
    figure(1),
    plot(task_vec(:,1), 'r', 'Linewidth', 1.5), hold on, plot(task_vec(:,10), 'r--', 'Linewidth', 1.5)
    legend('p_1', 'p ref_1', 'Location', 'best');
    title("Circular trajectory - EE cartesian X coordinate")
    xlabel('time [s]'),
    ylabel('EE position [m]'),
    grid on
    
    % EE position y-axis
    figure(2),
    plot(task_vec(:,2), 'b', 'Linewidth', 1.5), hold on, plot(task_vec(:,11), 'b--', 'Linewidth', 1.5)
    legend('p_2', 'p ref_2', 'Location', 'best');
    title("Circular trajectory - EE cartesian Y coordinate")
    xlabel('time [s]'),
    ylabel('EE position [m]'),
    grid on
    
    % EE position z-axis
    figure(3),
    plot(task_vec(:,3), 'g', 'Linewidth', 1.5), hold on, plot(task_vec(:,12), 'g--', 'Linewidth', 1.5)
    legend('p_2', 'p ref_2', 'Location', 'best');
    title("Circular trajectory - EE cartesian Z coordinate")
    xlabel('time [s]'),
    ylabel('EE position [m]'),
    grid on
end

%% cartesian errors

if show_error == 1
    
    % EE position error x-axis
    figure(4),
    plot(task_vec(:,1)-task_vec(:,10), 'r', 'Linewidth', 1.5)
    legend('e_x', 'Location', 'best');
    title("Circular trajectory - EE cartesian error - X coordinate")
    xlabel('time [s]'),
    ylabel('error [m]'),
    grid on
    
    % EE position error y-axis
    figure(5),
    plot(task_vec(:,2)-task_vec(:,11), 'b', 'Linewidth', 1.5)
    legend('e_y','Location', 'best');
    title("Circular trajectory - EE cartesian error - Y coordinate")
    xlabel('time [s]'),
    ylabel('error [m]'),
    grid on
    
    % EE position error z-axis
    figure(6),
    plot(task_vec(:,3)-task_vec(:,12), 'g', 'Linewidth', 1.5)
    legend('e_z', 'Location', 'best');
    title("Circular trajectory - EE cartesian error - Z coordinate")
    xlabel('time [s]'),
    ylabel('error [m]'),
    grid on
end

%% variance plot

if show_std == 1
    
    % prediction standard deviations
    figure(7),
    plot(sqrt(variances(:,1)), 'Linewidth', 1), hold on
    plot(sqrt(variances(:,2)), 'Linewidth', 1), hold on
    plot(sqrt(variances(:,3)), 'Linewidth', 1), hold on
    plot(sqrt(variances(:,4)), 'Linewidth', 1), hold on
    plot(sqrt(variances(:,5)), 'Linewidth', 1), hold on
    plot(sqrt(variances(:,6)), 'Linewidth', 1), hold on
    plot(sqrt(variances(:,7)), 'Linewidth', 1)
    
    legend('Std_1', 'Std_2', 'Std_3', 'Std_4', 'Std_5', 'Std_6', 'Std_7', 'Location', 'best');
    title("Circular trajectory - Prediction standard deviations")
    xlabel('time [s]'),
    ylabel('std'),
    grid on
    
end

%% scatter plot

if show_scatter == 1
    
    % scatter_plot, real vs desired
    figure(),
    sc1 = scatter3(task_vec(:,1), task_vec(:,2), task_vec(:,3));
    hold on
    sc2 = scatter3(task_vec(:,10), task_vec(:,11), task_vec(:,12));
    legend('Robot motion','Desired trajectory')
    set(sc1, 'Linewidth', 0.75)
    set(sc2, 'Linewidth', 0.75)
    title('KUKA EE Cartesian Actual vs Desired Path')
end

%% show robot

if show_animation == 1
    
    path = p_ref';
    
    % robot animation
    figure();
    show(kuka);
    
    view(3)
    ax = gca;
    ax.Projection = 'orthographic';
    hold on
    plot3(path(:,1), path(:,2), path(:,3), 'r','Linewidth',2)
    
    for i = 1 : 5 : length(t)-1
        show(kuka, q(i,:),'PreservePlot',false);
        plot3(task_vec(i,1), task_vec(i,2), task_vec(i,3), 'o', 'Color', 'b', 'Linewidth', 2)
        title(['Time = ' num2str(t(i))])
        drawnow
    end
    drawnow
end

%% save data
if save_data == 1
    model = 'model3';
    
    if gp_flag == 1
        if opt_flag == 1
            data_path = ['results/' model '/' traj '/linGP/opt' ];
        elseif opt_flag == 0
            data_path = ['results/' model '/' traj '/linGP/no_opt' ];
        end
        dlmwrite(fullfile(data_path, 'std.txt'), sqrt(variances))
        dlmwrite(fullfile(data_path, 'q.txt'), q)
        dlmwrite(fullfile(data_path, 'dq.txt'), dq)
        dlmwrite(fullfile(data_path, 'ddq.txt'), ddq)
        dlmwrite(fullfile(data_path, 'p.txt'), task_vec(:,1:3))
        dlmwrite(fullfile(data_path, 'dp.txt'), task_vec(:,4:6))
        dlmwrite(fullfile(data_path, 'ddp.txt'), task_vec(:,7:9))
        dlmwrite(fullfile(data_path, 'p_ref.txt'), task_vec(:,10:12))
        dlmwrite(fullfile(data_path, 'dp_ref.txt'), task_vec(:,13:15))
        dlmwrite(fullfile(data_path, 'ddp_ref.txt'), task_vec(:,16:18))
        dlmwrite(fullfile(data_path, 'corrective_torques.txt'), corrective_torques)
        dlmwrite(fullfile(data_path, 'needed_torques.txt'), needed_torques)
        dlmwrite(fullfile(data_path, 'singular_values.txt'), singular_values)
        dlmwrite(fullfile(data_path, 'torques.txt'), torques)
        disp('Data successfully saved!')
    elseif gp_flag == 0
        data_path = ['results/' model '/' traj '/no_learning/' ];
        dlmwrite(fullfile(data_path, 'q.txt'), q)
        dlmwrite(fullfile(data_path, 'dq.txt'), dq)
        dlmwrite(fullfile(data_path, 'ddq.txt'), ddq)
        dlmwrite(fullfile(data_path, 'p.txt'), task_vec(:,1:3))
        dlmwrite(fullfile(data_path, 'dp.txt'), task_vec(:,4:6))
        dlmwrite(fullfile(data_path, 'ddp.txt'), task_vec(:,7:9))
        dlmwrite(fullfile(data_path, 'p_ref.txt'), task_vec(:,10:12))
        dlmwrite(fullfile(data_path, 'dp_ref.txt'), task_vec(:,13:15))
        dlmwrite(fullfile(data_path, 'ddp_ref.txt'), task_vec(:,16:18))
        dlmwrite(fullfile(data_path, 'singular_values.txt'), singular_values)
        dlmwrite(fullfile(data_path, 'torques.txt'), torques)
        disp('Data successfully saved!')
    end
end

%% save workspace

if save_workspace == 1
    workspace_path = ['workspaces/' traj '/'];
    if gp_flag == 1
        if opt_flag == 1
            save([workspace_path 'linGP'])
        elseif opt_flag == 0
            save([workspace_path 'linGP_noOpt'])
        end
    elseif gp_flag == 0
        save([workspace_path 'no_learning'])
    end
    disp('Workspace saved!')
end