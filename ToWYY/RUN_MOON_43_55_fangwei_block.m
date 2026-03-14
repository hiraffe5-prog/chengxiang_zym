%%
close all;
clear;
clc;
load sub_region_dem38.mat
%% 参数设置
addpath(genpath('./子函数'));

Tr = 2.5e-3;                % 发射脉冲时宽
PRT = 16.666667e-3;         % 慢时间间隔
BW = 0.5e6;                 % 带宽
f0 = 3.12e9;              	% 雷达工作频率
c = 299792458;           	% 光速
lamda = c/f0;         	    % 波长
Tcoh=300;                   % 合成孔径时间
fs = 1e6;                   % 采样率
Nr = fs*Tr;                 % chirp信号采样点数
Nrg = round(PRT*fs);        % 雷达接收总采样点数
Kr = BW/Tr;                 % 调频斜率
Rmoon = 1737.40e3;          % 月球半径

Para.fs = fs;
Para.Nrg = Nrg;
Para.PRT = PRT;
Para.f0 = f0;
Para.Nr = Nr;
Para.Kr = Kr;
Para.Tr = Tr;
Para.c = c;
Para.lamda = lamda;
Para.BW = BW;
Para.Tcoh = Tcoh;
Para.Rmoon = Rmoon;
Para.B = ceil(11*60/Para.PRT);
%% 场景仿真

lat = 39;             % 目标中心点的纬度
lon = 51.5;             % 目标中心点的经度
lat_max = 1.5;         	% 整个场景沿纬度方向的最大范围。
lon_max = 1.5;         	% 整个场景沿经度方向的最大范围。

target_lat = linspace(lat-lat_max/2,lat+lat_max/2,384);     % 纬度向散射体
target_lon = linspace(lon-lon_max/2,lon+lon_max/2,384);     % 经度向散射体

num_target_lat = length(target_lat);        % 纬度向散射体点数
num_target_lon = length(target_lon);        % 经度向散射体点数


target_H = double(sub_region);      % 用来存放每个散射体的高度信息，矩阵。

lat_mat = target_lat.'*ones(1,num_target_lon);              %纬度网格化
lon_mat = ones(num_target_lat,1)*target_lon;                %经度网格化
target_x = (Rmoon+target_H).*cosd(lat_mat).*cosd(lon_mat);  %散射体的x坐标
target_y = (Rmoon+target_H).*cosd(lat_mat).*sind(lon_mat);  %散射体的y坐标
target_z = (Rmoon+target_H).*sind(lat_mat);                 %散射体的z坐标
center_x = (Rmoon).*cosd(lat).*cosd(lon);  %散射体的x坐标
center_y = (Rmoon).*cosd(lat).*sind(lon);  %散射体的y坐标
center_z = (Rmoon).*sind(lat);  
% 作图
figure;
[X,Y] = meshgrid(target_lon,target_lat);
surf(X,Y,target_H);
title('月球DEM');
clear X;clear Y;

Para.lat = lat;
Para.lon = lon;
Para.lat_max = lat_max;
Para.lon_max = lon_max;

folder = './echo_data8';

% 如果文件夹不存在，则创建它
if ~exist(folder, 'dir')
    mkdir(folder);
end
%% 星历读入
Path_Name = ['./数据/'];                                      %相对路径
File_Name = ['WGC_StateVector_20250629.csv'];          %雷达的MCMF坐标
FIle_Name_Full = [Path_Name File_Name];                      %全路径

start_time1 = '2025-06-29T01-28-54-000Z';       % 输入时间段
end_time1   = '2025-06-29T01-44-53-000Z';
[~,Radar_MCMF_xyz1] = nasa_crop_time(FIle_Name_Full,start_time1,end_time1);
x_Radar1 = Radar_MCMF_xyz1(:,1)*1000;     % 雷达在MCMF下的x坐标，m
y_Radar1 = Radar_MCMF_xyz1(:,2)*1000;     % 雷达在MCMF下的y坐标，m
z_Radar1 = Radar_MCMF_xyz1(:,3)*1000;     % 雷达在MCMF下的z坐标，m
Ttotal = length(x_Radar1);               % 插值前的方位向点数
Naz = ceil(Ttotal/PRT);                 % 插值后的方位向点数
ta = (0:1:Naz-1)*PRT;                   % 方位向时间
x_par_radar1 = polyfit(0:1:Ttotal-1,x_Radar1,3);  % 使用三次函数来对雷达轨迹进行拟合
y_par_radar1 = polyfit(0:1:Ttotal-1,y_Radar1,3);
z_par_radar1 = polyfit(0:1:Ttotal-1,z_Radar1,3);
x_radar1 = polyval(x_par_radar1,ta)';             % 根据拟合的三次方程进行插值
y_radar1 = polyval(y_par_radar1,ta)';
z_radar1 = polyval(z_par_radar1,ta)';
xyz_radar1 = [x_radar1 y_radar1 z_radar1];

start_time2 = '2025-06-29T02-06-54-000Z';      % 输入时间段
end_time2   = '2025-06-29T02-22-53-000Z';
[~,Radar_MCMF_xyz2] = nasa_crop_time(FIle_Name_Full,start_time2,end_time2);
x_Radar2 = Radar_MCMF_xyz2(:,1)*1000;     % 雷达在MCMF下的x坐标，m
y_Radar2 = Radar_MCMF_xyz2(:,2)*1000;     % 雷达在MCMF下的y坐标，m
z_Radar2 = Radar_MCMF_xyz2(:,3)*1000;     % 雷达在MCMF下的z坐标，m
x_par_radar2 = polyfit(0:1:Ttotal-1,x_Radar2,3);  % 使用三次函数来对雷达轨迹进行拟合
y_par_radar2 = polyfit(0:1:Ttotal-1,y_Radar2,3);
z_par_radar2 = polyfit(0:1:Ttotal-1,z_Radar2,3);
x_radar2 = polyval(x_par_radar2,ta)';             % 根据拟合的三次方程进行插值
y_radar2 = polyval(y_par_radar2,ta)';
z_radar2 = polyval(z_par_radar2,ta)';
xyz_radar2 = [x_radar2 y_radar2 z_radar2];
%% 分辨率计算
Trans = [-sind(lon) cosd(lon) 0;-sind(lat)*cosd(lon) -sind(lat)*sind(lon) cosd(lat);cosd(lat)*cosd(lon) cosd(lat)*sind(lon) sind(lat)];
vector1 = Trans*[x_radar1(1)-center_x;y_radar1(1)-center_y;z_radar1(1)-center_z];
vector2 = Trans*[x_radar1(end)-center_x;y_radar1(end)-center_y;z_radar1(end)-center_z];
vector1 = [vector1(1),vector1(2),0];
vector2 = [vector2(1),vector2(2),0];
dot_product = dot(vector1, vector2);
magnitude1 = norm(vector1);
magnitude2 = norm(vector2);
cos_theta = dot_product / (magnitude1 * magnitude2);
theta = acos(cos_theta);
rou_a1 = lamda/2/theta;
vector1 = Trans*[x_radar2(1)-center_x;y_radar2(1)-center_y;z_radar2(1)-center_z];
vector2 = Trans*[x_radar2(end)-center_x;y_radar2(end)-center_y;z_radar2(end)-center_z];
vector1 = [vector1(1),vector1(2),0];
vector2 = [vector2(1),vector2(2),0];
dot_product = dot(vector1, vector2);
magnitude1 = norm(vector1);
magnitude2 = norm(vector2);
cos_theta = dot_product / (magnitude1 * magnitude2);
theta = acos(cos_theta);
rou_a2 = lamda/2/theta;
disp(['第一轨的方位向分辨率', num2str(rou_a1)]);
disp(['第二轨的方位向分辨率', num2str(rou_a2)]);

vector1 = [x_radar1(end/2)-center_x,y_radar1(end/2)-center_y,z_radar1(end/2)-center_z];
vector2 = [center_x,center_y,center_z];
dot_product = dot(vector1, vector2);
norm1 = norm(vector1);
norm2 = norm(vector2);
theta_rad = acos(dot_product / (norm1 * norm2));
rou_d1 = c/2/BW/sin(theta_rad);
vector1 = [x_radar2(end/2)-center_x,y_radar2(end/2)-center_y,z_radar2(end/2)-center_z];
vector2 = [center_x,center_y,center_z];
dot_product = dot(vector1, vector2);
norm1 = norm(vector1);
norm2 = norm(vector2);
theta_rad = acos(dot_product / (norm1 * norm2));
rou_d2 = c/2/BW/sin(theta_rad);
disp(['第一轨的距离向分辨率', num2str(rou_d1)]);
disp(['第二轨的距离向分辨率', num2str(rou_d2)]);
%% 模糊数计算
R0 = 2*sqrt((x_radar1(end/2)-target_x(1,1))^2 + (y_radar1(end/2)-target_y(1,1))^2 + (z_radar1(end/2)-target_z(1,1))^2);    % 雷达某一时刻到场景某点的斜距，用于计算模糊数
Para.amb_num = floor(R0/Para.c/Para.PRT);           % 模糊数计算
Para.tr = (0:1:(Para.Nrg-1))*Para.PRT/Para.Nrg;     % 快时间轴，0-PRT，采样率为fs=Nrg/PRT
tr_mtx = ones(Naz,1)*Para.tr;    % 距离时间轴矩阵，大小：Naz*Nrg
Para.R0=R0;
Para.Nrg = Nrg;
Para.Naz = Naz;
Para.x_radar1=x_radar1;
Para.y_radar1=y_radar1;
Para.z_radar1=z_radar1;
Para.x_radar2=x_radar2;
Para.y_radar2=y_radar2;
Para.z_radar2=z_radar2;
%% DEM法向量逐点确定
[N, M] = size(target_x); % 获取矩阵的尺寸
% 初始化法向量矩阵
normals_x = zeros(N, M);  % 存储法向量的 x 分量
normals_y = zeros(N, M);  % 存储法向量的 y 分量
normals_z = zeros(N, M);  % 存储法向量的 z 分量
% 对于每个点计算法向量
for i = 2:N-1
    for j = 2:M-1
        E = [target_x(i-1,j),target_y(i-1,j),target_z(i-1,j)];  % 左边一点
        F = [target_x(i+1,j),target_y(i+1,j),target_z(i+1,j)];  % 右边一点
        G = [target_x(i,j-1),target_y(i,j-1),target_z(i,j-1)];  % 上边一点
        H = [target_x(i,j+1),target_y(i,j+1),target_z(i,j+1)];  % 下边一点
        FE = F-E;   % 左右向量
        HG = H-G;   % 上下向量
        normal = cross(FE, HG);         % 叉乘计算法向量
        normal_mag = norm(normal);      % 法向量模值
        normal = normal / normal_mag;  % 单位化法向量
        
        % 存储法向量的 x, y, z 分量
        normals_x(i, j) = normal(1);
        normals_y(i, j) = normal(2);
        normals_z(i, j) = normal(3);
    end
end

% 处理边缘点（边缘点使用单边差分法）
% 对于第一行、第一列、最后一行和最后一列进行单边差分法处理
for i = 1:N
    for j = 1:M
        if i == 1 || i == N || j == 1 || j == M
            % 处理边缘点的法向量，可以使用简单的单边差分计算
            if i == 1
                dx1 = target_x(i+1, j) - target_x(i, j);
                dy1 = target_y(i+1, j) - target_y(i, j);
                dz1 = target_z(i+1, j) - target_z(i, j);
                normal = [dx1, dy1, dz1];
            elseif i == N
                dx1 = target_x(i-1, j) - target_x(i, j);
                dy1 = target_y(i-1, j) - target_y(i, j);
                dz1 = target_z(i-1, j) - target_z(i, j);
                normal = [dx1, dy1, dz1];
            elseif j == 1
                dx1 = target_x(i, j+1) - target_x(i, j);
                dy1 = target_y(i, j+1) - target_y(i, j);
                dz1 = target_z(i, j+1) - target_z(i, j);
                normal = [dx1, dy1, dz1];
            elseif j == M
                dx1 = target_x(i, j-1) - target_x(i, j);
                dy1 = target_y(i, j-1) - target_y(i, j);
                dz1 = target_z(i, j-1) - target_z(i, j);
                normal = [dx1, dy1, dz1];
            end
            normal_mag = norm(normal);
            normal = normal / normal_mag;
            normals_x(i, j) = normal(1);
            normals_y(i, j) = normal(2);
            normals_z(i, j) = normal(3);
        end
    end
end
filename = fullfile(folder, 'echo_Para.mat');  % 拼接完整路径
save(filename, 'Para');  % 保存变量
%% 造回波
A_phase = 2*pi*rand(num_target_lat,num_target_lon);        % 随机相位
tic
tic
% 配置GPU
gpuDevice(1);
reset(gpuDevice);  % 确保GPU内存清空
% 转换常量到GPU变量
tr_mtx_gpu = gpuArray(tr_mtx);
f0_gpu = gpuArray(f0);
c_gpu = gpuArray(c);
Kr_gpu = gpuArray(Kr);
amb_num_gpu = gpuArray(Para.amb_num);
PRT_gpu = gpuArray(Para.PRT);
xyz_radarA = gpuArray(xyz_radar1(round(Naz/2),:));  % 预转换常用向量
xyz_radarB = gpuArray(xyz_radar2(round(Naz/2),:));

% 初始化GPU数组
one_Nrg_block = gpuArray.ones(1,Nrg);  % 用于后续分块计算

% 设置分块参数（根据显存调整）
block_size_a = round(60/PRT);
num_blocks_a = ceil(Naz / block_size_a);

for i = 1:num_blocks_a
    cols = (i-1)*block_size_a+1 : min(i*block_size_a, Naz);
    current_cols = cols(1):cols(end);
    s_echoA_local = gpuArray.zeros(length(current_cols), Nrg);
    s_echoB_local = gpuArray.zeros(length(current_cols), Nrg);
    filename = fullfile(folder,['echo' num2str(i) '.mat']);
    % 主循环优化
    for kk = 1:num_target_lat   % 遍历各个纬度向
        KK = kk;
        % 预转换目标坐标到CPU单精度
        target_x_kk = target_x(KK, :);  % 该纬度下散射体的x坐标
        target_y_kk = target_y(KK, :);  % 该纬度下散射体的y坐标
        target_z_kk = target_z(KK, :);  % 该纬度下散射体的z坐标
        normals_x_kk = normals_x(KK, :);
        normals_y_kk = normals_y(KK, :);
        normals_z_kk = normals_z(KK, :);
        A_phase_kk = A_phase(KK, :);
    
        for ll = 1:num_target_lon
            % 计算距离向量
            LL = ll;
            
            delta_tA = 2/c*sqrt((x_radar1(current_cols) - target_x_kk(LL)).^2 + ...
                    (y_radar1(current_cols) - target_y_kk(LL)).^2 + ...
                    (z_radar1(current_cols) - target_z_kk(LL)).^2);    %回波时延
            delta_tA_gpu = gpuArray(delta_tA);            %转化GPU变量
            delta_tB = 2/c*sqrt((x_radar2(current_cols) - target_x_kk(LL)).^2 + ...
                    (y_radar2(current_cols) - target_y_kk(LL)).^2 + ...
                    (z_radar2(current_cols) - target_z_kk(LL)).^2);    %回波时延    
            delta_tB_gpu = gpuArray(delta_tB);            %转化GPU变量    
            norm_temp = [normals_x_kk(LL), normals_y_kk(LL), normals_z_kk(LL)];
            norm_temp_gpu = gpuArray(norm_temp);
            cos_deltaA = abs(dot(norm_temp_gpu, xyz_radarA)) / ...
                    (vecnorm(norm_temp_gpu) * vecnorm(xyz_radarA));
            cos_deltaB = abs(dot(norm_temp_gpu, xyz_radarB)) / ...
                    (vecnorm(norm_temp_gpu) * vecnorm(xyz_radarB));
            Phase = A_phase_kk(LL);
            % 分块处理距离维度
            Delta_tA_block = delta_tA_gpu * one_Nrg_block;
            % 脉冲匹配条件
            time_diffA = tr_mtx_gpu(current_cols,:) + amb_num_gpu.*PRT_gpu - Delta_tA_block;% 接收信号快时间轴
            conditionA = abs(time_diffA) <= Tr/2; % 认为发射信号是-Tr/2-Tr/2，距离门
        
            % 相位计算（合并指数运算）
            phaseA = exp(1j * pi * (-2 * f0_gpu * Delta_tA_block + Kr_gpu * time_diffA.^2));
            % 当前块的信号贡献
            block_contributionA = conditionA .* phaseA;    % 门函数、chirp项与载频向（下变频后残留）
            % block_contributionA = conditionA .* phaseA  ;
            % 累加到总结果
            s_echoA_local = s_echoA_local + block_contributionA;
            
            Delta_tB_block = delta_tB_gpu * one_Nrg_block;
            time_diffB = tr_mtx_gpu(current_cols,:) + amb_num_gpu.*PRT_gpu - Delta_tB_block;% 接收信号快时间轴
            conditionB = abs(time_diffB) <= Tr/2; % 认为发射信号是-Tr/2-Tr/2，距离门
            phaseB = exp(1j * pi * (-2 * f0_gpu * Delta_tB_block + Kr_gpu * time_diffB.^2));
            block_contributionB = conditionB .* phaseB ;    % 门函数、chirp项与载频向（下变频后残留）
            % block_contributionB = conditionB .* phaseB ;
            s_echoB_local = s_echoB_local + block_contributionB;
        end
    end
    % 返回结果到CPU
    s_echoA = gather(s_echoA_local);
    s_echoB = gather(s_echoB_local); 
    save(filename,'s_echoA','s_echoB','-v7.3');
    % 每完成一个纬度释放CPU内存
    clear s_echoA s_echoB s_echoA_local s_echoB_local;
    disp(['已成像分钟数:' num2str(i) ]);
end
reset(gpuDevice);  % 清空GPU内存
toc
% colormap(gray)