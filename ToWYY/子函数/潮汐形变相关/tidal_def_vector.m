% 本函数的作用是计算MCMF坐标系下一定经纬度范围内的月球表面各点的潮汐形变向量
% 输入量
% 1.R：观测时的地月间距（km）
% 2.sep_lat, sep_lon: 观测时的地球星下点纬度与经度（°）
% 3.target_lat_range, target_lon_range: 需要分析潮汐形变的月面经纬度范围（°）
% 输出量
% 1.X_def, Y_def, Z_def: 潮汐形变在MCMF坐标系下的三维坐标矩阵（m）
% 2.target_lat, target_lon: 需要分析潮汐形变的月面经纬度矩阵（°）

function [X_def, Y_def, Z_def, target_lat, target_lon] = tidal_def_vector(R,sep_lat,sep_lon,target_lat_range,target_lon_range)
%% 计算相关参数设置
% 参数设置
mu_earth = 398600.435436;                   % 地球地心引力常数GM（km3/s2）
R_moon = 1737.4;                            % 月球半径（km）
g = 1.623e-3;                               % 月球表面重力常数（km/s2）
h2 = 0.0424;                                % 垂直位移勒夫数
l2 = 0.0107;                                % 水平位移勒夫数
% 参数计算
[target_lon,target_lat] = meshgrid(target_lon_range,target_lat_range);                  % 构建待分析区域经纬度网格（°）

%% 计算雷达星下点与各待分析位置在MCMF坐标系下的单位向量
% 地球星下点单位向量
ux_ = cosd(sep_lat)*cosd(sep_lon)*ones(size(target_lon)); uy_ = cosd(sep_lat)*sind(sep_lon)*ones(size(target_lon)); uz_ = sind(sep_lat)*ones(size(target_lon));
% 待分析位置单位向量
ux = cosd(target_lat).*cosd(target_lon); uy = cosd(target_lat).*sind(target_lon); uz = sind(target_lat);

%% 计算潮汐形变值
% 计算潮汐形变公式中各个部分参数
alpha_cos = ux*ux_ + uy*uy_ + uz*uz_;                                                                               % 对应为地球星下点单位向量与待分析位置单位向量之积                                                
a = ( ( mu_earth.*(R_moon.^2) ) ./ ( 2.*g.*(R.^3) ) ) .* (3*(alpha_cos.^2)-1) .* h2 .* 1000 .*100;                  % 潮汐形变中垂直形变的系数，单位为cm                                              
b = ( ( 3.*mu_earth.*(R_moon^2) ) ./ ( g.*(R.^3) ) ) .* alpha_cos .* l2 .* 1000 .*100;                              % 潮汐形变中水平形变的系数，单位为cm
c = alpha_cos;                                                                                                      % 水平形变方向组合系数
% 计算MCMF坐标系下潮汐形变的三维坐标矩阵
X_def = -a.*ux-b.*(ux_-c.*ux); Y_def = -a.*uy-b.*(uy_-c.*uy); Z_def = -a.*uz-b.*(uz_-c.*uz);                           % 单位为cm
X_def = X_def./100; Y_def = Y_def./100; Z_def = Z_def./100;                                                         % 单位为m
end