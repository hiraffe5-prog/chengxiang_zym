%%%%%%%%%%%%%%%%%%%%%% 基于相干系数的三级图像配准方法 %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%         Author：LG           %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%      Data：2020.09.1         %%%%%%%%%%%%%%%%%%%%
%————————————————————————————————————————————————%
%——————————————————————更新状态——————————————————————%
% 版 本 号：1.0.0
% 近期修改：无
% 执 行 人：李莫凡
% 审 核 人：无
%————————————————————————————————————————————————%
%——————————————————————函数介绍——————————————————————%
% 函数功能：将副图像信息配准到主图像，以实现后续的干涉处理
%————————————————————————————————————————————————%
%                    参数名             参数定义              数据尺寸          数据类型
% 输入参数：         master              主图像           二维矩阵（M*N）        double
%                    slave              副图像           二维矩阵（M*N）        double
%             sample_block_length     采样窗大小          单个元素（1*1）       double
%             search_block_length     搜索窗大小          单个元素（1*1）       double
%                threshold_cor       相干系数阈值         单个元素（1*1）       double
%————————————————————————————————————————————————%
%                    参数名             参数定义              数据尺寸          数据类型
% 输出参数：        slave_reg           配准后图像         二维矩阵（M*N）       double
%————————————————————————————————————————————————%
function [slave_reg]=Register_main_pure_v2(master,slave,search_row_col,sample_block_length,search_block_length,threshold_cor,overlapping,upsample_row_col)
%% 
% 第一步：校正图像整体像素偏移
% 第二步：校正图像像素级几何畸变
% 第三步：校正图像亚像素级几何积畸变
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Na,Nr] = size(master);    
%% 粗配准：估计图像整体偏移量(这里辅图像是只能往上移吗，如果是的话就是有问题的)
    sample_block_row=Na;
    sample_block_col=Nr;
    search_row=search_row_col;
    search_col=search_row_col;
    search_block_row=Na+search_row;
    search_block_col=Nr+search_col;
    master_temp=zeros(Na+search_row,Nr+search_col);
    slave_temp=zeros(Na+search_row,Nr+search_col);
%     master_temp(search_row+1:Na+search_row,search_col+1:Nr+search_col)=master;
%     slave_temp(search_row+1:Na+search_row,search_col+1:Nr+search_col)=slave;
    master_temp(search_row/2+1:Na+search_row/2,search_col/2+1:Nr+search_col/2)=master;
    slave_temp(search_row/2+1:Na+search_row/2,search_col/2+1:Nr+search_col/2)=slave;
    [delta_row,delta_col] = Phase_cor(master_temp,slave_temp,Na+search_row,Nr+search_col,sample_block_row,sample_block_col,search_block_row,search_block_col);
    save SaveData/Phase_cor.mat delta_row delta_col
    %%
    if delta_row>0
        slave_temp1=[slave(delta_row+1:Na,:);zeros(delta_row,Nr)]; % 向上平移delta_row
    else
        slave_temp1=[zeros(-delta_row,Nr);slave(1:Na+delta_row,:)];   % 向下平移delta_row
    end
    if delta_col>0
        slave_temp1=[slave_temp1(:,delta_col+1:Nr) zeros(Na,delta_col)]; % 向左平移delta_col
    else
        slave_temp1=[zeros(Na,-delta_col) slave_temp1(:,1:delta_col+Nr)];    % 向右平移delta_col
    end
   
 %% 亚像素级配准
    sample_block_row = sample_block_length;
    sample_block_col = sample_block_length;
    search_block_row = search_block_length;
    search_block_col = search_block_length; % 经验：图像尺寸/搜索窗尺寸在4.5到6比较好
    overlapping = overlapping;              % 控制点间隔1/overlapping个样本窗宽度
    upsample_row = upsample_row_col;
    upsample_col = upsample_row_col;
    %
    [row_pix_master,column_pix_master,row_offset,col_offset,max_cor] = Match_subpix(master,slave_temp1,Na,Nr,...
                       sample_block_row,sample_block_col,search_block_row,search_block_col,upsample_row,upsample_col,overlapping);
    save SaveData/Match_subpix.mat row_pix_master column_pix_master row_offset col_offset max_cor
    %}
    % 对配准像素进行二维拟合
    load SaveData/Match_subpix.mat
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