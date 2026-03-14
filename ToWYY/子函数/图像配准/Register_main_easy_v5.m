% 本函数的作用是利用之前获取的不同分辨率的配准参数进行图像配准

function [slave_reg]=Register_main_easy_v5(master,slave,Para2D_m,flag_pl,flag_spl)
[Na,Nr] = size(master);                                                                                         
slave_temp1 = slave;
% 加载配准参数
if flag_pl == 1
    load SaveData/Phase_cor_0.mat delta_row delta_col
    
else
    delta_row = 0; delta_col = 0;
end
if flag_spl == 1
    load SaveData/Match_subpix_0_add.mat row_offset_all col_offset_all
else
    row_offset_all = zeros(Na,Nr); col_offset_all = zeros(Na,Nr);
end
row_offset_all = row_offset_all+delta_row; col_offset_all = col_offset_all+delta_col;
% 修正偏移量，使之等比例扩大到当前的坐标系（因为我们知道两次网格的经纬度范围是相同的，点数都是等比例增加的）
[Na_,Nr_] = size(row_offset_all);               % 获得配准参数对应的Na和Nr
row_offset_all = row_offset_all*(Na/Na_);       % 偏移量修正
col_offset_all = col_offset_all*(Nr/Nr_);
% 构建配准参数与当前情况各点对应的经纬度值（经纬度轴）
lon_min = min(Para2D_m.lon); lon_max = max(Para2D_m.lon);
lat_min = min(Para2D_m.lat); lat_max = max(Para2D_m.lat);
lon_axis = linspace(lon_min,lon_max,Na); lat_axis = linspace(lat_min,lat_max,Nr);
[lon2D,lat2D] = meshgrid(lon_axis,lat_axis);
lon_axis_ = linspace(lon_min,lon_max,Na_); lat_axis_ = linspace(lat_min,lat_max,Nr_);
[lon2D_,lat2D_] = meshgrid(lon_axis_,lat_axis_);
% 将偏移量按照经纬度轴进行重采样
row_offset_all = interp2(lon2D_,lat2D_,row_offset_all,lon2D,lat2D);
col_offset_all = interp2(lon2D_,lat2D_,col_offset_all,lon2D,lat2D);
% 将偏移量值转化为配准后的新位置
xx = row_offset_all+(1:Na)'*ones(1,Nr);
yy = col_offset_all+ones(Na,1)*(1:Nr);
% 将原数据配准到指定区域并进行双线性插值
add=500;
slave_add=zeros(Na+add,Nr+add);
slave_add(add/2+1:Na+add/2,add/2+1:Nr+add/2)=slave_temp1;
% 双线性插值，权值矩阵
offset_matrix1=zeros(Na,Nr);
offset_matrix2=zeros(Na,Nr);
offset_matrix3=zeros(Na,Nr);
offset_matrix4=zeros(Na,Nr);
% 数据矩阵
matrix1=zeros(Na,Nr);
matrix2=zeros(Na,Nr);
matrix3=zeros(Na,Nr);
matrix4=zeros(Na,Nr);
for ii=1:Na
  %权值矩阵
  offset_matrix1(ii,:)=1-(xx(ii,:)-floor(xx(ii,:))); % 上权
  offset_matrix2(ii,:)=xx(ii,:)-floor(xx(ii,:));    %下权
  offset_matrix3(ii,:)=1-(yy(ii,:)-floor(yy(ii,:)));  %左权
  offset_matrix4(ii,:)=yy(ii,:)-floor(yy(ii,:));     %右权
  % 数据矩阵
  matrix1(ii,:)=slave_add(sub2ind(size(slave_add),floor(xx(ii,:)+add/2),floor(yy(ii,:)+add/2)));      % 左上数据
  matrix2(ii,:)=slave_add(sub2ind(size(slave_add),floor(xx(ii,:)+add/2),floor(yy(ii,:)+add/2+1)));    % 右上数据
  matrix3(ii,:)=slave_add(sub2ind(size(slave_add),floor(xx(ii,:)+add/2+1),floor(yy(ii,:)+add/2)));    % 左下数据
  matrix4(ii,:)=slave_add(sub2ind(size(slave_add),floor(xx(ii,:)+add/2+1),floor(yy(ii,:)+add/2+1)));  % 右下数据    
end
slave_reg=matrix1.*offset_matrix1.*offset_matrix3+matrix2.*offset_matrix1.*offset_matrix4+...
                matrix3.*offset_matrix2.*offset_matrix3+matrix4.*offset_matrix2.*offset_matrix4;
end