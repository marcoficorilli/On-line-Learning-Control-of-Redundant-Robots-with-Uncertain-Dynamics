%% This file returns and saves in the proper repository (that must exist!) several plots from different experiments for each task.

% Note: the save commands are all commented.

clear all
close all
clc

% utilities
model = 'model3';   % ONLY MODEL 3 AVAILABLE! - { model1, model2, model3 } -> depending on how much uncertainty we want to add on the masses
traj = 'heart';             % { circular, sinusoidal, heart}

% set up the save folder location
path_fig = ['figures/' traj '/fig/'];
path_eps = ['figures/' traj '/eps/'];

if (strcmp(traj, 'circular'))
    traj_name = 'Circular';
elseif (strcmp(traj, 'sinusoidal'))
    traj_name = 'Sinusoidal';
elseif (strcmp(traj, 'heart'))
    traj_name = 'Heart';
end

t = 0:10e-3:5;
t = t(2:end)';

% (auto) set path
data_path1 = ['results/' model '/' traj '/linGP/opt'];
data_path2 = ['results/' model '/' traj '/linGP/no_opt'];
data_path3 = ['results/' model '/' traj '/no_learning'];

% load data
p1 = load([data_path1 '/p.txt']);              % EE cartesian ACTUAL position
p2 = load([data_path2 '/p.txt']);
p3 = load([data_path3 '/p.txt']);
p_ref1 = load([data_path1 '/p_ref.txt']);      % EE cartesian REFERENCE position
p_ref2 = load([data_path2 '/p_ref.txt']); 
p_ref3 = load([data_path3 '/p_ref.txt']); 
torques1 = load([data_path1 '/torques.txt']);   % torques learning + optimization
torques2 = load([data_path2 '/torques.txt']);   % torques learning without optimization
torques3 = load([data_path3 '/torques.txt']);   % torques no learning

corrective_torques1 = load([data_path1 '/corrective_torques.txt']);
corrective_torques2 = load([data_path2 '/corrective_torques.txt']);
needed_torques = load([data_path1 '/needed_torques.txt']);

std1 = load([data_path1 '/std.txt']);
std2 = load([data_path2 '/std.txt']);

std1_norm = [];
std2_norm = [];
for i = 1:size(std1,1)
    std1_norm = vertcat(std1_norm, norm(std1(i,:),2));
    std2_norm = vertcat(std2_norm, norm(std2(i,:),2));
end

% set min and max value to scale torque plots
if min(min(torques1)) < min(min(torques2))
    tau_min1 = min(min(torques1));
else
    tau_min1 = min(min(torques2));
end

if max(max(torques1)) > max(max(torques2))
    tau_max1 = max(max(torques1));
else
    tau_max1 = max(max(torques2));
end

if min(min(torques1)) < min(min(torques3))
    tau_min2 = min(min(torques1));
else
    tau_min2 = min(min(torques3));
end

if max(max(torques1)) > max(max(torques3))
    tau_max2 = max(max(torques1));
else
    tau_max2 = max(max(torques3));
end

if tau_min1 < tau_min2
    tau_min = tau_min1 - 1;
else 
    tau_min = tau_min2 - 1;
end

if tau_max1 > tau_max2
    tau_max = tau_max1 + 1;
else
    tau_max = tau_max2 + 1;
end

if max(max(std1)) > max(max(std2))
    std_max = max(max(std1)) + 0.5;
else
    std_max = max(max(std2)) + 0.5;
end


%% error comparison opt vs not opt
f1 = figure();
f1.WindowState = 'maximized';
subplot(3,1,1)
plot(t, abs(p1(:,1)-p_ref1(:,1)), 'r', 'Linewidth', 1.5), hold on
plot(t, abs(p2(:,1)-p_ref2(:,1)), 'r--', 'Linewidth', 1.5)
legend('learning + optimization | e_x |', 'learning | e_x |', 'Location', 'best', 'FontSize', 8);
title("X coordinate error")
xlabel('time [s]'),
ylabel('error [m]'),
grid on

subplot(3,1,2)
plot(t, abs(p1(:,2)-p_ref1(:,2)), 'b', 'Linewidth', 1.5), hold on
plot(t, abs(p2(:,2)-p_ref2(:,2)), 'b--', 'Linewidth', 1.5)
legend('learning + optimization | e_y |', 'learning | e_y |', 'Location', 'best', 'FontSize', 8);
title("Y coordinate error")
xlabel('time [s]'),
ylabel('error [m]'),
grid on

subplot(3,1,3)
plot(t, abs(p1(:,3)-p_ref1(:,3)), 'g', 'Linewidth', 1.5), hold on
plot(t, abs(p2(:,3)-p_ref2(:,3)), 'g--', 'Linewidth', 1.5)
legend('learning + optimization | e_z |', 'learning | e_z |', 'Location', 'best', 'FontSize', 8);
title("Z coordinate error")
xlabel('time [s]'),
ylabel('error [m]'),
grid on

sgtitle([traj_name ' - Linearized GP - Optimized vs Non-Optimized'], 'Color', 'black', 'FontSize', 14);

% saveas(gcf, [path_fig 'error_opt_vs_noOpt'], 'fig')
% saveas(gcf, [path_eps 'error_opt_vs_noOpt'], 'epsc')
% 
% disp('error_opt_vs_noOpt saved!')
% disp(' ')

%% error comparison learning vs no learning
f2 = figure();
f2.WindowState = 'maximized';
subplot(3,1,1)
plot(t, abs(p1(:,1)-p_ref1(:,1)), 'r', 'Linewidth', 1.5), hold on
plot(t, abs(p3(:,1)-p_ref3(:,1)), 'r--', 'Linewidth', 1.5)
legend('learning | e_x |', 'model based | e_x |', 'Location', 'best', 'FontSize', 8);
title("X coordinate error")
xlabel('time [s]'),
ylabel('error [m]'),
grid on

subplot(3,1,2)
plot(t, abs(p1(:,2)-p_ref1(:,2)), 'b', 'Linewidth', 1.5), hold on
plot(t, abs(p3(:,2)-p_ref3(:,2)), 'b--', 'Linewidth', 1.5)
legend('learning | e_y |', 'model based | e_y |', 'Location', 'best', 'FontSize', 8);
title("Y coordinate error")
xlabel('time [s]'),
ylabel('error [m]'),
grid on

subplot(3,1,3)
plot(t, abs(p1(:,3)-p_ref1(:,3)), 'g', 'Linewidth', 1.5), hold on
plot(t, abs(p3(:,3)-p_ref3(:,3)), 'g--', 'Linewidth', 1.5)
legend('learning | e_z |', 'model based | e_z |', 'Location', 'best', 'FontSize', 8);
title("Z coordinate error")
xlabel('time [s]'),
ylabel('error [m]'),
grid on

sgtitle([traj_name ' - Linearized GP - Learning vs Pure model based'], 'Color', 'black', 'FontSize', 14);

% saveas(gcf, [path_fig 'error_learning_vs_noLearning'], 'fig')
% saveas(gcf, [path_eps 'error_learning_vs_noLearning'], 'epsc')
% 
% disp('error_learning_vs_noLearning saved!')
% disp(' ')


%% Torque comparison
figure()
plot(t, torques1(:,1), 'Linewidth', 1), hold on
plot(t, torques1(:,2), 'Linewidth', 1), hold on
plot(t, torques1(:,3), 'Linewidth', 1), hold on
plot(t, torques1(:,4), 'Linewidth', 1), hold on
plot(t, torques1(:,5), 'Linewidth', 1), hold on
plot(t, torques1(:,6), 'Linewidth', 1), hold on
plot(t, torques1(:,7), 'Linewidth', 1)
legend('Tau_1', 'Tau_2', 'Tau_3', 'Tau_4', 'Tau_5', 'Tau_6', 'Tau_7', 'Location', 'best', 'FontSize', 8, 'Orientation', 'horizontal', 'NumColumns', 3);
ylim([tau_min tau_max])
title([traj_name ' - Commanded Torques - Learning + Optimization'])
xlabel('time [s]'),
ylabel('tau [Nm]'),
grid on
% saveas(gcf, [path_fig 'commanded_torques_opt'], 'fig')
% saveas(gcf, [path_eps 'commanded_torques_opt'], 'epsc')
% disp('torques_opt saved!')
% disp(' ')


figure()
plot(t, torques2(:,1), 'Linewidth', 1), hold on
plot(t, torques2(:,2), 'Linewidth', 1), hold on
plot(t, torques2(:,3), 'Linewidth', 1), hold on
plot(t, torques2(:,4), 'Linewidth', 1), hold on
plot(t, torques2(:,5), 'Linewidth', 1), hold on
plot(t, torques2(:,6), 'Linewidth', 1), hold on
plot(t, torques2(:,7), 'Linewidth', 1)
legend('Tau_1', 'Tau_2', 'Tau_3', 'Tau_4', 'Tau_5', 'Tau_6', 'Tau_7', 'Location', 'best', 'FontSize', 8, 'Orientation', 'horizontal', 'NumColumns', 3);
ylim([tau_min tau_max])
title([traj_name ' - Commanded Torques - Learning' ])
xlabel('time [s]'),
ylabel('tau [Nm]'),
grid on
% saveas(gcf, [path_fig 'commanded_torques_noOpt'], 'fig')
% saveas(gcf, [path_eps 'commanded_torques_noOpt'], 'epsc')
% disp('torques_noOpt saved!')
% disp(' ')

figure()
plot(t, torques3(:,1), 'Linewidth', 1), hold on
plot(t, torques3(:,2), 'Linewidth', 1), hold on
plot(t, torques3(:,3), 'Linewidth', 1), hold on
plot(t, torques3(:,4), 'Linewidth', 1), hold on
plot(t, torques3(:,5), 'Linewidth', 1), hold on
plot(t, torques3(:,6), 'Linewidth', 1), hold on
plot(t, torques3(:,7), 'Linewidth', 1)
legend('Tau_1', 'Tau_2', 'Tau_3', 'Tau_4', 'Tau_5', 'Tau_6', 'Tau_7', 'Location', 'best', 'FontSize', 8, 'Orientation', 'horizontal', 'NumColumns', 3);
ylim([tau_min tau_max])
title([traj_name ' - Commanded Torques - Model Based' ])
xlabel('time [s]'),
ylabel('tau [Nm]'),
grid on
% saveas(gcf, [path_fig 'commanded_torques_noLearning'], 'fig')
% saveas(gcf, [path_eps 'commanded_torques_noLearning'], 'epsc')
% disp('torques_noLearning saved!')
% disp(' ')

%% std comparison
figure()
plot(t, std1(:,1), 'Linewidth', 1), hold on
plot(t, std1(:,2), 'Linewidth', 1), hold on
plot(t, std1(:,3), 'Linewidth', 1), hold on
plot(t, std1(:,4), 'Linewidth', 1), hold on
plot(t, std1(:,5), 'Linewidth', 1), hold on
plot(t, std1(:,6), 'Linewidth', 1), hold on
plot(t, std1(:,7), 'Linewidth', 1)
legend('Std_1', 'Std_2', 'Std_3', 'Std_4', 'Std_5', 'Std_6', 'Std_7', 'Location', 'best', 'FontSize', 8);
ylim([0 std_max])
title([traj_name ' - Prediction Standard Deviation - Optimized'])
xlabel('time [s]'),
ylabel('std  [rad/s^2]'),
grid on
% saveas(gcf, [path_fig 'std_opt'], 'fig')
% saveas(gcf, [path_eps 'std_opt'], 'epsc')
% disp('std_opt saved!')
% disp(' ')

figure()
plot(t, std2(:,1), 'Linewidth', 1), hold on
plot(t, std2(:,2), 'Linewidth', 1), hold on
plot(t, std2(:,3), 'Linewidth', 1), hold on
plot(t, std2(:,4), 'Linewidth', 1), hold on
plot(t, std2(:,5), 'Linewidth', 1), hold on
plot(t, std2(:,6), 'Linewidth', 1), hold on
plot(t, std2(:,7), 'Linewidth', 1)
legend('Std_1', 'Std_2', 'Std_3', 'Std_4', 'Std_5', 'Std_6', 'Std_7', 'Location', 'best', 'FontSize', 8);
ylim([0 std_max])
title([traj_name ' - Prediction Standard Deviation - Not Optimized'])
xlabel('time [s]'),
ylabel('std  [rad/s^2]'),
grid on
% saveas(gcf, [path_fig 'std_noOpt'], 'fig')
% saveas(gcf, [path_eps 'std_noOpt'], 'epsc')
% disp('std_noOpt saved!')
% disp(' ')
%%
f3 = figure();
f3.WindowState = 'maximized';
plot(t, std1_norm, 'r', 'Linewidth', 1.5), hold on,
plot(t, std2_norm, 'r--', 'Linewidth', 1.5), 
legend('Optimized', 'Non Optimized', 'Location', 'best', 'FontSize', 8);
title([traj_name ' - Prediction Standard Deviation - Norm Comparison'])
xlabel('time [s]'),
ylabel('std  [rad/s^2]'),
grid on
% saveas(gcf, [path_fig 'std_norm'], 'fig')
% saveas(gcf, [path_eps 'std_norm'], 'epsc')
% disp('std_norm saved!')
% disp(' ')

%%
f4 = figure();
f4.WindowState = 'maximized';
subplot(4,1,1)
plot(t, corrective_torques1(:,1), 'b', 'Linewidth', 1.5), hold on
plot(t, corrective_torques2(:,1), 'r', 'Linewidth', 1.5)
plot(t, needed_torques(:,1), 'g', 'Linewidth', 1.5)
legend('learning + optimization', 'learning', 'needed torques', 'Location', 'best', 'FontSize', 12, 'Orientation', 'horizontal');
title('1^{st} joint')
xlabel('time [s]'),
ylabel('torque  [Nm]'),
grid on

subplot(4,1,2)
plot(t, corrective_torques1(:,2), 'b', 'Linewidth', 1.5), hold on
plot(t, corrective_torques2(:,2), 'r', 'Linewidth', 1.5)
plot(t, needed_torques(:,2), 'g', 'Linewidth', 1.5)
title('2^{nd} joint')
xlabel('time [s]'),
ylabel('torque  [Nm]'),
grid on

subplot(4,1,3)
plot(t, corrective_torques1(:,3), 'b', 'Linewidth', 1.5), hold on
plot(t, corrective_torques2(:,3), 'r', 'Linewidth', 1.5)
plot(t, needed_torques(:,3), 'g', 'Linewidth', 1.5)
title('3^{rd} joint')
xlabel('time [s]'),
ylabel('torque  [Nm]'),
grid on

subplot(4,1,4)
plot(t, corrective_torques1(:,4), 'b', 'Linewidth', 1.5), hold on
plot(t, corrective_torques2(:,4), 'r', 'Linewidth', 1.5)
plot(t, needed_torques(:,4), 'g', 'Linewidth', 1.5)
title('4^{th} joint')
xlabel('time [s]'),
ylabel('torque  [Nm]'),
grid on
sgtitle([traj_name ' - Corrective Torques Comparison - Joints 1-4'], 'Color', 'black', 'FontSize', 14);
% saveas(gcf, [path_fig 'corr_torques_1-4'], 'fig')
% saveas(gcf, [path_eps 'corr_torques_1-4'], 'epsc')
% 
% disp('corr_torques_1-4 saved!')
% disp(' ')
%%
%
%
%
f5 = figure();
f5.WindowState = 'maximized';
subplot(3,1,1)
plot(t, corrective_torques1(:,5), 'b', 'Linewidth', 1.5), hold on
plot(t, corrective_torques2(:,5), 'r', 'Linewidth', 1.5)
plot(t, needed_torques(:,5), 'g', 'Linewidth', 1.5)
title('5^{th} joint')
xlabel('time [s]'),
ylabel('torque  [Nm]'),
legend('learning + optimization', 'learning', 'needed torques', 'Location', 'best', 'FontSize', 12, 'Orientation', 'horizontal');
grid on

subplot(3,1,2)
plot(t, corrective_torques1(:,6), 'b', 'Linewidth', 1.5), hold on
plot(t, corrective_torques2(:,6), 'r', 'Linewidth', 1.5)
plot(t, needed_torques(:,6), 'g', 'Linewidth', 1.5)
title('6^{th} joint')
xlabel('time [s]'),
ylabel('torque  [Nm]'),
grid on

subplot(3,1,3)
plot(t, corrective_torques1(:,7), 'b', 'Linewidth', 1.5), hold on
plot(t, corrective_torques2(:,7), 'r', 'Linewidth', 1.5)
plot(t, needed_torques(:,7), 'g', 'Linewidth', 1.5)
title('7^{th} joint')
xlabel('time [s]'),
ylabel('torque  [Nm]'),
grid on

sgtitle([traj_name ' - Corrective Torques Comparison - Joints 5-7'], 'Color', 'black', 'FontSize', 14);

% saveas(gcf, [path_fig 'corr_torques_5-7'], 'fig')
% saveas(gcf, [path_eps 'corr_torques_5-7'], 'epsc')
% 
% disp('corr_torques_5-7 saved!')
% disp(' ')

%%
% scatter_plot, real vs desired

f7 = figure();
f7.WindowState = 'maximized';
sc1 = scatter3(p1(:,1), p1(:,2), p1(:,3), [], 'blue');
hold on
sc2 = scatter3(p3(:,1), p3(:,2), p3(:,3), [], 'red');
hold on
sc3 = scatter3(p_ref1(:,1), p_ref1(:,2), p_ref1(:,3), [], 'green');
set(sc1, 'Linewidth', 1.75)
set(sc2, 'Linewidth', 1.75)
set(sc3, 'Linewidth', 1.75)
title('Comparison Cartesian Actual vs Desired Path')
legend('Learning', 'Model Based', 'Desired')

% saveas(gcf, [path_fig 'scatter'], 'fig')
% saveas(gcf, [path_eps 'scatter'], 'epsc')
% 
% disp('scatter saved!')
% disp(' ')
