clear,clc;
addpath(genpath('./子函数'));

%% 1.参数设置
load echo_MOON_43_55_points0.mat % 回波参数
load sub_region_dem_2.mat       % 场景高程
% load echo_MOON_43_55_single_tar.mat
target_H = double(sub_region);  % 读取场景高程
Para.B = ceil(60/Para.PRT);     % 基线长度
Para.s_echoA = s_echoA;         % 回波

Para.lat_net = linspace(Para.lat-Para.lat_max,Para.lat+Para.lat_max,512);   % 成像网格纬度向
Para.lon_net = linspace(Para.lon-Para.lon_max,Para.lon+Para.lon_max,512);   % 成像网格经度向

lat_mat = Para.lat_net.'*ones(1,length(Para.lon_net));      %网格化
lon_mat = ones(length(Para.lat_net),1)*Para.lon_net;

target_lat = linspace(Para.lat-Para.lat_max/2,Para.lat+Para.lat_max/2,512);     % 纬度向散射体
target_lon = linspace(Para.lon-Para.lon_max/2,Para.lon+Para.lon_max/2,512);     % 经度向散射体
num_target_lat = length(target_lat);        % 纬度向散射体点数
num_target_lon = length(target_lon);        % 经度向散射体点数
lat_mat1 = target_lat.'*ones(1,num_target_lon);              %纬度网格化
lon_mat1 = ones(num_target_lat,1)*target_lon;                %经度网格化
Para.tar_x = (Para.Rmoon+target_H).*cosd(lat_mat1).*cosd(lon_mat1);   % 含地形网格的x坐标
Para.tar_y = (Para.Rmoon+target_H).*cosd(lat_mat1).*sind(lon_mat1);   % 含地形网格的y坐标
Para.tar_z = (Para.Rmoon+target_H).*sind(lat_mat1);                  % 含地形网格的z坐标

    
Para.flat_x = (Para.Rmoon-3000).*cosd(lat_mat).*cosd(lon_mat);           % 不含地形网格的x坐标
Para.flat_y = (Para.Rmoon-3000).*cosd(lat_mat).*sind(lon_mat);           % 不含地形网格的y坐标
Para.flat_z = (Para.Rmoon-3000).*sind(lat_mat);                          % 不含地形网格的z坐标

flag_reg_pixel_level = 1;                                               % 决定本次处理是否进行粗配准的标记，取1为进行取0为不进行
search_row_col = 100;                                                    % 搜索窗相比原图像扩展的行数列数
% 设置精配准参数
flag_reg_sub_pixel_level = 1;                                           % 决定本次处理是否进行精配准的标记，取1为进行取0为不进行
% sample_window = 64;                                                    % 样本窗大小
% search_window = 80;                                                    % 搜索窗大小
sample_window = 48;                                                    % 样本窗大小
search_window = 128;                                                    % 搜索窗大小
hist_num = 30;                                                          % 偏移量统计分析时直方图个数
overlapping = 4;                                                        % 样本窗重叠倍数
upsample_row_col = 20;                                                  % 行列向升采样倍数
% 设置相位滤波参数
filter_window = 21;                                                      % 滤波窗口大小
% 设置经纬度再截取参数
flag_image_cut = 0;                                             
%% 1.bp成像
disp('1.开始bp成像')
[s_imag_A] = BPA_imaging(Para,1);
[s_imag_B] = BPA_imaging(Para,2); 
disp('bp成像完成')
%% 2.图像配准与干涉图生成
% 图像配准模块待补充
load imag_A_points.mat
s_imag_A = s_ac;
load imag_B_points.mat
s_imag_B = s_ac;
figure;
imagesc(abs(s_imag_A));
xlabel('经度');
ylabel('纬度');
figure;
imagesc(abs(s_imag_B));
xlabel('经度');
ylabel('纬度');
% lon_range = [42.5,44.5];
% lat_range = [54.5,56.5];
% lon_lat_range=[lon_range,lat_range];
% RD_geocoding3(s_imag_A,round((Para.Naz-Para.B)/2),lon_lat_range,del_lon,del_lat,Para.pos_xyz,Para.Tr,Para.Naz,f0);
% RD_geocoding3(s_imag_B,round((Para.Naz+Para.B)/2),lon_lat_range,del_lon,del_lat,Para.pos_xyz,Para.Tr,Para.Naz,f0);
disp('2.开始图像配准与干涉')
% [s_imag_B_col] = Register_main_pure_v4(s_imag_A,s_imag_B,search_row_col,sample_window,...
%         search_window,hist_num,overlapping,upsample_row_col,flag_reg_pixel_level,...
%         flag_reg_sub_pixel_level,0);     % 进行图像配准，双线性插值（多视后图像的配准结果在存储时使用标号0）
s_imag = s_imag_A.*conj(s_imag_B);
figure;imagesc(angle(s_imag));title('保含地形相位');
disp('图像配准与干涉图生成完成')

%% 3.平地相位去除
disp('3.开始进行平地相位去除')
PHY_flat_earth = calculate_Phase_flat(Para);
s_PHY_flat_earth = exp(1j*PHY_flat_earth);
s_after_flat_earth = s_imag.*conj(s_PHY_flat_earth);% 去平地相位后的干涉图（包括幅度和相位）
figure;imagesc(angle(s_PHY_flat_earth));title('理论计算得到的平地相位');
figure;imagesc(Para.lon_net,Para.lat_net,angle(s_after_flat_earth));title('去平地相位后的相位图');
xlabel('经度');
ylabel('纬度');

disp('平地与地形相位去除完成')

%% 4.相位滤波、相位解缠与高程反演

disp('4.开始进行形变反演')
PHY_after_filtering = Average_Filtering(angle(s_after_flat_earth),5,5);
figure;imagesc(Para.lon_net,Para.lat_net,PHY_after_filtering);title('滤波后的相位图');
disp('形变反演完成')

rmpath(genpath('./子函数'));