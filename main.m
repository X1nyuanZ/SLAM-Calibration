clear;clc;close all;

%% Main
land_marks = [-1 -1 -1;
              4 -4 -1;
              -0.5 2.5 0;
              4.5 3.5 0;];
n_l = size(land_marks,1);
dt=0.1;
t_span=[0:dt:15];
Q_in_Euler = [pi/6;0;pi/3;0;0;1];
R_c = eul2rotm(Q_in_Euler(1:3,:)',"XYZ");
t_c = Q_in_Euler(4:6,:);
Q = [R_c t_c;0 0 0 1];

w = [0;0;-0.4];
v = [0;1;0];


gauss_noise=wgn(3*n_l,25000,-30);
noise = reshape(gauss_noise,[3,n_l,25000]);
% noise = zeros([3,n_l,25000]);
measurement_rate = 0.1;
[T_1,T_2,hat_T_1,hat_T_2,hat_z] = slam_calibration(land_marks,t_span,Q_in_Euler, ...
    noise,measurement_rate,w,v);

%% Data Prcocessing

data_length = size(hat_T_1,3);
[hat_p_1,hat_p_2,p_1,p_2,p_2_in_i] = deal(zeros(data_length,3));

hat_p_1_in_i=zeros(data_length,3);
hat_p_2_in_L1=zeros(data_length,4);
hat_p_2_in_i=zeros(data_length,4);
hat_z_2_in_L1=zeros(data_length,4*n_l);
hat_z_2_in_i=zeros(data_length,4*n_l);
hat_z_1_in_i=zeros(data_length,3*n_l);

hat_T_2_in_L1=zeros([4,4,data_length]);
hat_Q_reproduce=zeros([4,4,data_length]);

for i=1:1:data_length
    hat_p_1(i,:)=hat_T_1(1:3,4,i)';
    hat_p_2(i,:)=hat_T_2(1:3,4,i)';
    p_1(i,:)=T_1(1:3,4,i)';
    p_2(i,:)=T_2(1:3,4,i)';

    p_2_in_i(i,:) = (R_c*p_2(i,:)'+t_c)';
    hat_p_2_in_L1(i,:) = (hat_Q(:,:,i)*[hat_p_2(i,:)';1])';
    hat_z_2_in_L1(i,1:4*n_l)=reshape(hat_Q(:,:,i)*[reshape(hat_z(i,3*n_l+1:2*3*n_l),3,n_l);ones(1,n_l)],1,4*n_l);

    hat_p_2_in_i(i,:) = (T_inv(wave_T_1)*hat_p_2_in_L1(i,:)')';
end

for i=1:1:data_length
%     wave_T_1 = wave_T_his(1:4,1:4,i);
    R_1_residual=wave_T_1(1:3,1:3);
    t_1_residual=wave_T_1(1:3,4);
    hat_p_1_in_i(i,:)=(R_1_residual'*(hat_p_1(i,:)'-t_1_residual))';
    hat_z_1_in_i(i,1:3*n_l)=reshape((R_1_residual'*(reshape(hat_z(i,1:3*n_l),3,n_l)-t_1_residual)),1,3*n_l);
    hat_z_2_in_i(i,1:4*n_l)=reshape(T_inv(wave_T_1)*reshape(hat_z_2_in_L1(i,1:4*n_l),4,n_l),1,4*n_l);

    hat_T_2_in_L1(:,:,i) = hat_Q(:,:,i)*hat_T_2(:,:,i);
    hat_Q_reproduce(:,:,i) = T_inv(hat_T_1(:,:,i))*hat_T_2_in_L1(:,:,i);
end

wave_T_1_inv = T_inv(wave_T_1);

start_point = [0 0 0;0 0 1];
end_point = [hat_p_1_in_i(end,:);hat_p_2_in_i(end,1:3)];

hat_Q_reproduce_euler = rotm2eul(hat_Q_reproduce(1:3,1:3,:),'XYZ');
hat_Q_reproduce_position = reshape(hat_Q_reproduce(1:3,4,:),3,[])';

p_1_in_L1 = (wave_T_1*[p_1 ones(data_length,1)]')';
p_2_in_L1 = (wave_T_1*[p_2_in_i(:,1:3) ones(data_length,1)]')';

start_point_L1 = [p_1_in_L1(1,1:3);p_2_in_L1(1,1:3)];
end_point_L1 = [hat_p_1(end,:);hat_p_2_in_L1(1:3,end)'];
landmark_end_L1 = [hat_z(end,1:3);
    hat_z(end,4:6);
    hat_z(end,7:9);
    hat_z(end,10:12)];

for i=1:1:data_length
    norm_e_Q(i,1) = norm(Q-hat_Q_reproduce(:,:,i));
end
%% Fusion Map（LiDAR-1 Frame）

figure
hold on
grid on
view([-25 25])
xlabel('$x$','FontSize',15,'Interpreter','latex');
ylabel('$y$','FontSize',15,'Interpreter','latex');
zlabel('$z$','FontSize',15,'Interpreter','latex');

p1 = plot3(hat_p_1(:,1),hat_p_1(:,2),hat_p_1(:,3),'b.-');
p2 = plot3(hat_p_2_in_L1(:,1)',hat_p_2_in_L1(:,2)',hat_p_2_in_L1(:,3)','k.-');

for i=1:1:n_l
    p3 = plot3(hat_z(:,3*(i-1)+1),hat_z(:,3*(i-1)+2),hat_z(:,3*(i-1)+3),'m-');
    p4 = plot3(hat_z_2_in_L1(:,4*(i-1)+1),hat_z_2_in_L1(:,4*(i-1)+2),hat_z_2_in_L1(:,4*(i-1)+3),'g-');
end

plot3(p_1_in_L1(:,1),p_1_in_L1(:,2),p_1_in_L1(:,3),'r');
plot3(p_2_in_L1(:,1),p_2_in_L1(:,2),p_2_in_L1(:,3),'r');
s1 = scatter3(landmark_end_L1(:,1),landmark_end_L1(:,2),landmark_end_L1(:,3),'^', ...
        'MarkerEdgeColor','k','MarkerFaceColor','m');
s2 = scatter3(start_point_L1(:,1),start_point_L1(:,2),start_point_L1(:,3),'p', ...
        'MarkerEdgeColor','k','MarkerFaceColor','r');
s3 = scatter3(end_point_L1(:,1),end_point_L1(:,2),end_point_L1(:,3),'o', ...
        'MarkerEdgeColor','k','MarkerFaceColor','r');
s4 = scatter3(0,0,0,'s', ...
    'MarkerEdgeColor','k','MarkerFaceColor','b');

title("LiDAR-1 Frame")

lg1 = legend([p1 p2 p3 p4],"$\hat{p}_{1}$","$\hat{Q}\bar{\hat{p}}_{2}$", ...
    "$\hat{z}_{j}^{1}$","$\hat{Q}\bar{\hat{z}}_{j}^{2}$", ...
    'Location','northeast','FontSize',12, ...
    'FontName','Cambria Math', 'Interpreter','latex');
ah1=axes('position',get(gca,'position'),'visible','off');
lg2 = legend(ah1,[s1 s2 s3 s4],"Landmark","Start point", ...
    "End point","Initial estimate", ...
    'Location','south','FontSize',12, ...
    'FontName','Cambria Math', 'Interpreter','latex');
set(lg2,'position',[0.8 0.4 0.1 0.1])

hold off

%% Fusion Map（World Frame）

figure
hold on
grid on
xlim([-1.5 6])
ylim([-4.5 4])
zlim([-1 1.5])
xlabel('$x$ (m)','FontSize',15,'Interpreter','latex');
ylabel('$y$ (m)','FontSize',15,'Interpreter','latex');
zlabel('$z$ (m)','FontSize',15,'Interpreter','latex');
view([-55 25])

p1 = plot3(hat_p_1_in_i(:,1),hat_p_1_in_i(:,2),hat_p_1_in_i(:,3),'b.-');
p2 = plot3(hat_p_2_in_i(:,1)',hat_p_2_in_i(:,2)',hat_p_2_in_i(:,3)','k.-');

for i=1:1:n_l
    p3 = plot3(hat_z_1_in_i(:,3*(i-1)+1),hat_z_1_in_i(:,3*(i-1)+2),hat_z_1_in_i(:,3*(i-1)+3),'m-');
    p4 = plot3(hat_z_2_in_i(:,4*(i-1)+1),hat_z_2_in_i(:,4*(i-1)+2),hat_z_2_in_i(:,4*(i-1)+3),'g-');
end

plot3(p_1(:,1),p_1(:,2),p_1(:,3),'r')
plot3(p_2_in_i(:,1),p_2_in_i(:,2),p_2_in_i(:,3),'r')

s1 = scatter3(land_marks(:,1),land_marks(:,2),land_marks(:,3),'^', ...
        'MarkerEdgeColor','k','MarkerFaceColor','m');
s2 = scatter3(start_point(:,1),start_point(:,2),start_point(:,3),'p', ...
        'MarkerEdgeColor','k','MarkerFaceColor','r');
s3 = scatter3(end_point(:,1),end_point(:,2),end_point(:,3),'o', ...
        'MarkerEdgeColor','k','MarkerFaceColor','r');
s4 = scatter3(wave_T_1_inv(1,4),wave_T_1_inv(2,4),wave_T_1_inv(3,4),'s', ...
        'MarkerEdgeColor','k','MarkerFaceColor','b');

title("World Frame")

lg1 = legend([p1 p2 p3 p4],"$\hat{p}_{1}$","$\hat{Q}\bar{\hat{p}}_{2}$", ...
    "$\hat{z}_{j}^{1}$","$\hat{Q}\bar{\hat{z}}_{j}^{2}$", ...
    'Location','northeast','FontSize',12, ...
    'FontName','Cambria Math', 'Interpreter','latex');
ah1=axes('position',get(gca,'position'),'visible','off');
lg2 = legend(ah1,[s1 s2 s3 s4],"Landmark","Start point", ...
    "End point","Initial estimate", ...
    'Location','south','FontSize',12, ...
    'FontName','Cambria Math', 'Interpreter','latex');
set(lg2,'position',[0.8 0.4 0.1 0.1])

hold off

%% Extrinsic Parmeter \hat{Q} (Error, Rotation, Translation)

figure

subplot(3,1,1)
plot((1:1:data_length)*dt,norm_e_Q,'b.-','LineWidth',1.5,'MarkerSize',10);
grid on
xlim([0,15])
ylim([0,2.5])
ylabel('$\|Q-\hat{Q}\|$','FontSize',15,'Interpreter','latex');
title("Error")

subplot(3,1,2)
plot((1:1:data_length)*dt,hat_Q_reproduce_euler,'.-','LineWidth',1.5,'MarkerSize',10);
grid on
xlim([0,15])
% xlabel('time','FontSize',15,'Interpreter','latex');
ylabel('Euler (rad)','FontSize',13,'Interpreter','latex');
title("Rotation")

subplot(3,1,3)
plot((1:1:data_length)*dt,hat_Q_reproduce_position,'.-','LineWidth',1.5,'MarkerSize',10);
grid on
xlim([0,15])
ylim([-1,2])
xlabel('time (s)','FontSize',13,'Interpreter','latex');
ylabel('Translation (m)','FontSize',13,'Interpreter','latex');
title("Translation")

function out_matrix = T_inv(in_matrix)
    M_1 = in_matrix(1:3,1:3);
    m_2 = in_matrix(1:3,4);
    out_matrix = [M_1' -M_1'*m_2;0 0 0 1];
end