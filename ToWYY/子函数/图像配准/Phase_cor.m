%%%%%%%%%%%%%%%%%%%%%% 基于相干系数的三级图像配准方法 %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%        配准参数子函数         %%%%%%%%%%%%%%%%%%%%
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
% 函数功能：计算粗配准所需的配准参数，即图像整体偏移量信息
%————————————————————————————————————————————————%
%                    参数名             参数定义              数据尺寸          数据类型
% 输入参数：       data_master           主图像           二维矩阵（M*N）        double
%                  data_slave           副图像           二维矩阵（M*N）        double
%                     Na/Nr         图像尺寸（列/行）      单个元素（1*1）       double
%             sample_block_row/col     采样窗大小         单个元素（1*1）        double
%             search_block_row/col     搜索窗大小         单个元素（1*1）        double
%————————————————————————————————————————————————%
%                    参数名             参数定义              数据尺寸          数据类型
% 输出参数：       delta_row       整体横坐标偏移量       单个元素（1*1）       double
%                 delta_col       整体纵坐标偏移量       单个元素（1*1）       double
%————————————————————————————————————————————————%
function [delta_row,delta_col] = Phase_cor(data_master,data_slave,Na,Nr,...
          sample_block_row,sample_block_col,search_block_row,search_block_col)
%UNTITLED Summary of this function goes here

%%
%% 数据及参数读取
master_image_new = data_master;             % 这里的主辅图像都预先在外围扩展了一圈搜索窗与样本窗边长差距的一半，便于搜索窗的构建
slave_image_new = data_slave;
row_new = Na;
col_new = Nr;
clear data_master;
clear data_slave;
%% 确定窗口个数
% 根据紧凑样本窗个数确定整体窗口个数
% 样本窗口与搜索窗口个数相同（与某个样本窗对应的搜索窗是比其大一圈的一个窗口，代表预先设定的样本窗可能平移的范围）
% 注意：边缘部分必须留出每个样本窗对应搜索窗的搜索区域
block_number_row = floor((row_new-(search_block_row-sample_block_row))/sample_block_row); 
block_number_col = floor((col_new-(search_block_col-sample_block_col))/sample_block_col); 
%% sample窗口控制点（左上角）坐标
%行列起始坐标相同（就是样本窗的左上角）
sample_control_pix_start_row = search_block_row/2-sample_block_row/2+1; 
sample_control_pix_start_col = search_block_col/2-sample_block_col/2+1; 
%行列尾部坐标不同（因原始数据行列数的限制）（就是样本窗的右下角）
sample_control_pix_end_row = block_number_row*sample_block_row+search_block_row/2-sample_block_row/2;
sample_control_pix_end_col = block_number_col*sample_block_col+search_block_col/2-sample_block_col/2;
% sample控制点向量（所有样本窗左上角的行号列号组成的向量）
line_sample_control_point_row = sample_control_pix_start_row:sample_block_row:sample_control_pix_end_row;
line_sample_control_point_col = sample_control_pix_start_col:sample_block_col:sample_control_pix_end_col;
% sample控制点矩阵(旧)
sample_control_point_row_old = repmat(line_sample_control_point_row',1,block_number_col);
sample_control_point_col_old = repmat(line_sample_control_point_col,block_number_row,1);
%%  sample样本点坐标矩阵(旧)
% 根据sample左上角点的坐标确定其右下方sample方格中心点的坐标，
% 此坐标输出后只用于二维拟合，所以不必要全部为整数
  row_pix_master=sample_control_point_row_old+(sample_block_row-1)/2;
  column_pix_master=sample_control_point_col_old+(sample_block_col-1)/2;
%% search窗口控制点（左上角）坐标
% 对应于每一个样本窗，搜索窗控制点一定在其左边（或上边）,其差值一定是窗口边长差的1/2
% search控制点向量
line_search_control_point_row = line_sample_control_point_row-(search_block_row/2-sample_block_row/2);
line_search_control_point_col = line_sample_control_point_col-(search_block_col/2-sample_block_col/2);

%% 
% 偏移量矩阵
row_offset = zeros(block_number_row,block_number_col);
col_offset  = zeros(block_number_row,block_number_col);
% 样本窗控制点再搜索窗中的相对坐标，依据此坐标为基准，求取偏移量
start_row_in_matrix = search_block_row/2-sample_block_row/2+1;
start_col_in_matrix = search_block_col/2-sample_block_col/2+1;

%% 滑窗
% 滑窗范围及坐标
% *因样本窗控制点为左上角点，偏移量最大值为窗口差值
% *使得求偏移量的起始参考点为左上角点
% *所以此处的坐标最小值固定为(1,1)
% *后续在截图确定峰值点时不用考虑参考坐标补偿问题
max_offset_window_row = floor((search_block_row-sample_block_row)/2);       % 配准最大可能进行的行向偏移
max_offset_window_col = floor((search_block_col-sample_block_col)/2);
max_row_in_matrix = start_row_in_matrix + max_offset_window_row;
max_col_in_matrix = start_col_in_matrix + max_offset_window_col;
%主图像0矩阵模板
master_board = zeros(search_block_row,search_block_col);
% 0及1构成的模板矩阵
ones_matirx = zeros(search_block_row,search_block_col);
ones_matirx(1:sample_block_row,1:sample_block_col)=1;       % 这个01模板里是在搜索窗左上角取的1
% 记录最大相关系数的值矩阵
max_cor = zeros(block_number_row,block_number_col);
%%

for mm = 1 : block_number_row
    for nn = 1 : block_number_col
         % 显示运行百分比
%          ((mm-1)*block_number_col+nn)/(block_number_row*block_number_col)*100
         % 控制点坐标
         m1 = line_sample_control_point_row(mm); % master左上角行坐标
         m2 = line_sample_control_point_col(nn);  % master左上角列坐标
         s1 = line_search_control_point_row(mm);   % slave左上角行坐标
         s2 = line_search_control_point_col(nn);    % slave左上角列坐标    
        
        % 取值
        master_block = master_image_new( m1:m1+sample_block_row-1 , m2 : m2+sample_block_col-1);    % 主图像取样本窗
        slave_block = slave_image_new( s1:s1+search_block_row-1 , s2 : s2+search_block_col-1);      % 辅图像取搜索窗
          % 实相干系数分子 直接主从图像频域滑窗 
        % 实相干系数的fft等于主从图像的fft之后再共轭相乘
        master_board(1:sample_block_row,1:sample_block_col) = master_block;
        Temp1 = abs(ifft2(conj(fft2(abs(master_board))).*fft2(abs(slave_block))));
        % 直接计算
        Temp2 = sqrt(sum(sum(abs(master_block).^2))).*...
                        sqrt(abs(ifft2(conj(fft2(abs(ones_matirx))).*fft2(abs(slave_block).^2))));
        % 实相干系数
        temp_cor = Temp1./Temp2;
    
        %求取偏移量(限制滑窗范围)
        part_temp_cor = temp_cor(1:max_row_in_matrix,1:max_col_in_matrix);
        
        %注意截图后的坐标起始值问题
        [row_new_in_matrix,col_new_in_matrix] = find(part_temp_cor==max(max(part_temp_cor)));
        row_offset(mm,nn) = row_new_in_matrix(1) - start_row_in_matrix;
        col_offset(mm,nn) = col_new_in_matrix(1) - start_col_in_matrix ;
        
         %记录最大相关系数的值
         max_cor(mm,nn) = max(max(part_temp_cor));           
    end
end
delta_row=row_offset;
delta_col=col_offset;
end
