% 本函数的作用是将成像结果与同经纬度范围内的DEM进行配准
function [DEM_reg] = Register_main_DEM(image,DEM,search_row_col,sample_block_length,search_block_length,threshold_cor,overlapping,upsample_row_col)
%% 参数获取
[Na,Nr] = size(image);
%% 参数预处理
% 对SAR图像数据，转化为幅度图，像素点值归一化
image = mat2gray(20*log10(abs(image)));
% 对DEM数据，求纬度方向梯度，像素点值归一化
DEM_ = zeros(Na,Nr);
DEM_(1:(end-1),:) = DEM(1:(end-1),:)-DEM(2:end,:);
DEM_ = mat2gray(DEM_);
% 为了增加特征对于两组数据的支配力度，将两者像素值低于平均值的值置为0
image(find(image>mean(mean(image))))=0;
DEM_(find(DEM_>mean(mean(DEM_))))=0;

%% 像素级配准
% 移动参数计算
sample_block_row=Na;
sample_block_col=Nr;
search_row=search_row_col;
search_col=search_row_col;
search_block_row=Na+search_row;
search_block_col=Nr+search_col;
master_temp=zeros(Na+search_row,Nr+search_col);
slave_temp=zeros(Na+search_row,Nr+search_col);
master_temp(search_row/2+1:Na+search_row/2,search_col/2+1:Nr+search_col/2)=image;
slave_temp(search_row/2+1:Na+search_row/2,search_col/2+1:Nr+search_col/2)=DEM_;
[delta_row,delta_col] = Phase_cor(master_temp,slave_temp,Na+search_row,Nr+search_col,sample_block_row,sample_block_col,search_block_row,search_block_col);
save(['./SaveData/Phase_cor_DEM.mat'],'delta_row','delta_col')
% 进行图像平移
if delta_row > 0
    DEM_reg = [DEM(delta_row+1:Na,:);zeros(delta_row,Nr)];
    DEM_ = [DEM_(delta_row+1:Na,:);zeros(delta_row,Nr)];
else
    DEM_reg = [zeros(-delta_row,Nr);DEM(1:Na+delta_row,:)];
    DEM_ = [zeros(-delta_row,Nr);DEM_(1:Na+delta_row,:)];
end
if delta_col > 0
    DEM_reg = [DEM_reg(:,delta_col+1:Nr) zeros(Na,delta_col)];
    DEM_ = [DEM_(:,delta_col+1:Nr) zeros(Na,delta_col)];
else
    DEM_reg = [zeros(Na,-delta_col) DEM_reg(:,1:delta_col+Nr)];
    DEM_ = [zeros(Na,-delta_col) DEM_(:,1:delta_col+Nr)];
end

 %% 亚像素级配准
sample_block_row = sample_block_length;
sample_block_col = sample_block_length;
search_block_row = search_block_length;
search_block_col = search_block_length; % 经验：图像尺寸/搜索窗尺寸在4.5到6比较好
overlapping = overlapping;              % 控制点间隔1/overlapping个样本窗宽度
upsample_row = upsample_row_col;
upsample_col = upsample_row_col;
if exist('SaveData/Match_subpix_DEM.mat','file') == 0
    %
    [row_pix_master,column_pix_master,row_offset,col_offset,max_cor] = Match_subpix(image,DEM_,Na,Nr,...
                       sample_block_row,sample_block_col,search_block_row,search_block_col,upsample_row,upsample_col,overlapping);
    save SaveData/Match_subpix.mat row_pix_master column_pix_master row_offset col_offset max_cor
    %}
else
    disp('完成DEM与SAR图像精配准时使用了上次保存的配准参数！')
    load SaveData/Match_subpix_DEM.mat
end
% 对配准像素进行二维拟合
m=size(row_pix_master);
max_cor_line=reshape(max_cor,m(1)*m(2),1);
L_old=reshape(row_pix_master,m(1)*m(2),1);
P_old=reshape(column_pix_master,m(1)*m(2),1);
L_new=reshape((row_pix_master+row_offset),m(1)*m(2),1);
P_new=reshape((column_pix_master+col_offset),m(1)*m(2),1);

mm=find(max_cor_line<threshold_cor);
L_old(mm)=[];
P_old(mm)=[];
L_new(mm)=[];
P_new(mm)=[];
% 拟合
x=L_old;
y=P_old;
z_row=L_new;
%拟合方程
fitobjectx=fit([x,y],z_row,'poly22');
z_column=P_new;
%拟合方程
fitobjecty=fit([x,y],z_column,'poly22');
% 对图像进行重采样
xx=zeros(Na,Nr);
yy=zeros(Na,Nr);
for ii=1:Na
    xx(ii,:)=fitobjectx(ii,1:Nr);
    yy(ii,:)=fitobjecty(ii,1:Nr);
end
% 将原数据配准到指定区域并进行双线性插值
add=2000;
slave_add=zeros(Na+add,Nr+add);
slave_add(add/2+1:Na+add/2,add/2+1:Nr+add/2)=DEM_reg;
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
DEM_reg=matrix1.*offset_matrix1.*offset_matrix3+matrix2.*offset_matrix1.*offset_matrix4+...
                matrix3.*offset_matrix2.*offset_matrix3+matrix4.*offset_matrix2.*offset_matrix4;
end