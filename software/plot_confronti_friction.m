%% This file returns and saves in the proper repository (that must exist!) all plots referred to the friction case for the sinusoidal task.

% Note: the save commands are all commented.

clear all
close all
clc

% utilities
model = 'model3';  
traj = 'sinusoidal';             % { sinusoidal}

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
data_path1 = ['results/' model '/' traj '/linGP/opt/friction'];
data_path2 = ['results/' model '/' traj '/no_learning/friction'];

% load data
p1 = load([data_path1 '/p.txt']);              % EE cartesian ACTUAL position
p2 = load([data_path2 '/p.txt']);
p_ref1 = load([data_path1 '/p_ref.txt']);      % EE cartesian REFERENCE position
p_ref2 = load([data_path2 '/p_ref.txt']); 

corrective_torques = load([data_path1 '/corrective_torques.txt']);
needed_torques_masses = load([data_path1 '/needed_torques_masses.txt']);
delta_eta = load([data_path1 '/delta_eta.txt']);

std1 = load([data_path1 '/std.txt']);

%% error comparison learning vs no learning
f1 = figure();
f1.WindowState = 'maximized';
subplot(3,1,1)
plot(t, abs(p1(:,1)-p_ref1(:,1)), 'r', 'Linewidth', 1.5), hold on
plot(t, abs(p2(:,1)-p_ref2(:,1)), 'r--', 'Linewidth', 1.5)
legend('learning | e_x |', 'model based | e_x |', 'Location', 'best', 'FontSize', 8);
title("X coordinate error")
xlabel('time [s]'),
ylabel('error [m]'),
grid on

subplot(3,1,2)
plot(t, abs(p1(:,2)-p_ref1(:,2)), 'b', 'Linewidth', 1.5), hold on
plot(t, abs(p2(:,2)-p_ref2(:,2)), 'b--', 'Linewidth', 1.5)
legend('learning | e_y |', 'model based | e_y |', 'Location', 'best', 'FontSize', 8);
title("Y coordinate error")
xlabel('time [s]'),
ylabel('error [m]'),
grid on

subplot(3,1,3)
plot(t, abs(p1(:,3)-p_ref1(:,3)), 'g', 'Linewidth', 1.5), hold on
plot(t, abs(p2(:,3)-p_ref2(:,3)), 'g--', 'Linewidth', 1.5)
legend('learning | e_z |', 'model based | e_z |', 'Location', 'best', 'FontSize', 8);
title("Z coordinate error")
xlabel('time [s]'),
ylabel('error [m]'),
grid on

sgtitle([traj_name ' - Linearized GP - Learning vs Pure model based'], 'Color', 'black', 'FontSize', 14);

% saveas(gcf, [path_fig 'error_learning_vs_noLearning_friction'], 'fig')
% saveas(gcf, [path_eps 'error_learning_vs_noLearning_friction'], 'epsc')
% 
% disp('error_learning_vs_noLearning_friction saved!')
% disp(' ')



%% std comparison
% figure()
% plot(t, std1(:,1), 'Linewidth', 1), hold on
% plot(t, std1(:,2), 'Linewidth', 1), hold on
% plot(t, std1(:,3), 'Linewidth', 1), hold on
% plot(t, std1(:,4), 'Linewidth', 1), hold on
% plot(t, std1(:,5), 'Linewidth', 1), hold on
% plot(t, std1(:,6), 'Linewidth', 1), hold on
% plot(t, std1(:,7), 'Linewidth', 1)
% legend('Std_1', 'Std_2', 'Std_3', 'Std_4', 'Std_5', 'Std_6', 'Std_7', 'Location', 'best', 'FontSize', 8);
% ylim([0 std_max])
% title([traj_name ' - Prediction Standard Deviation - Optimized'])
% xlabel('time [s]'),
% ylabel('std  [rad/s^2]'),
% grid on
% saveas(gcf, [path_fig 'std_opt'], 'fig')
% saveas(gcf, [path_eps 'std_opt'], 'epsc')
% disp('std_opt saved!')
% disp(' ')

%%
t2 = t(2:end-1);

f4 = figure();
f4.WindowState = 'maximized';
subplot(4,1,1)
plot(t2, corrective_torques(2+1:end,1) - needed_torques_masses(2+1:end,1), 'b', 'Linewidth', 1.5), hold on
plot(t2, delta_eta(2+1:end,1), 'r', 'Linewidth', 1.5)
legend('estimated', 'real', 'Location', 'best', 'FontSize', 12, 'Orientation', 'horizontal');
title('1^{st} joint')
xlabel('time [s]'),
ylabel('torque  [Nm]'),
grid on

subplot(4,1,2)
plot(t2, corrective_torques(2+1:end,2) - needed_torques_masses(2+1:end,2), 'b', 'Linewidth', 1.5), hold on
plot(t2, delta_eta(2+1:end,2), 'r', 'Linewidth', 1.5)
title('2^{nd} joint')
xlabel('time [s]'),
ylabel('torque  [Nm]'),
grid on

subplot(4,1,3)
plot(t2, corrective_torques(2+1:end,3) - needed_torques_masses(2+1:end,3), 'b', 'Linewidth', 1.5), hold on
plot(t2, delta_eta(2+1:end,3), 'r', 'Linewidth', 1.5)
title('3^{rd} joint')
xlabel('time [s]'),
ylabel('torque  [Nm]'),
grid on

subplot(4,1,4)
plot(t2,corrective_torques(2+1:end,4) - needed_torques_masses(2+1:end,4), 'b', 'Linewidth', 1.5), hold on
plot(t2, delta_eta(2+1:end,4), 'r', 'Linewidth', 1.5)
title('4^{th} joint')
xlabel('time [s]'),
ylabel('torque  [Nm]'),
grid on
sgtitle(['Sinusoidal - \Delta\eta estimate - Joints 1-4'], 'Color', 'black', 'FontSize', 14);

% saveas(gcf, [path_fig 'friction_joints_1_4'], 'fig')
% saveas(gcf, [path_eps 'friction_joints_1_4'], 'epsc')
% 
% disp('friction_joints_1_4 saved!')
% disp(' ')

f5 = figure();
f5.WindowState = 'maximized';
subplot(3,1,1)
plot(t2, corrective_torques(2+1:end,5) - needed_torques_masses(2+1:end,5), 'b', 'Linewidth', 1.5), hold on
plot(t2, delta_eta(2+1:end,5), 'r', 'Linewidth', 1.5)
title('5^{th} joint')
xlabel('time [s]'),
ylabel('torque  [Nm]'),
legend('estimated', 'real', 'Location', 'best', 'FontSize', 12, 'Orientation', 'horizontal');
grid on

subplot(3,1,2)
plot(t2, corrective_torques(2+1:end,6) - needed_torques_masses(2+1:end,6), 'b', 'Linewidth', 1.5), hold on
plot(t2, delta_eta(2+1:end,6), 'r', 'Linewidth', 1.5)
title('6^{th} joint')
xlabel('time [s]'),
ylabel('torque  [Nm]'),
grid on

subplot(3,1,3)
plot(t2, corrective_torques(2+1:end,7) - needed_torques_masses(2+1:end,7) , 'b', 'Linewidth', 1.5), hold on
plot(t2, delta_eta(2+1:end,7), 'r', 'Linewidth', 1.5)
title('7^{th} joint')
xlabel('time [s]'),
ylabel('torque  [Nm]'),
grid on

sgtitle(['Sinusoidal - \Delta\eta estimate - Joints 5-7'], 'Color', 'black', 'FontSize', 14);
% saveas(gcf, [path_fig 'friction_joints_5_7'], 'fig')
% saveas(gcf, [path_eps 'friction_joints_5_7'], 'epsc')
% 
% disp('friction_joints_5_7 saved!')
% disp(' ')

%%
% scatter_plot, real vs desired

f7 = figure();
f7.WindowState = 'maximized';
sc1 = scatter3(p1(:,1), p1(:,2), p1(:,3), [], 'blue');
hold on
sc2 = scatter3(p2(:,1), p2(:,2), p2(:,3), [], 'red');
hold on
sc3 = scatter3(p_ref1(:,1), p_ref1(:,2), p_ref1(:,3), [], 'green');
set(sc1, 'Linewidth', 1.75)
set(sc2, 'Linewidth', 1.75)
set(sc3, 'Linewidth', 1.75)
title('Comparison Cartesian Actual vs Desired Path')
legend('Learning', 'Model Based', 'Desired')

% saveas(gcf, [path_fig 'friction_scatter'], 'fig')
% saveas(gcf, [path_eps 'friction_scatter'], 'epsc')
% 
% disp('friction_scatter saved!')
% disp(' ')
