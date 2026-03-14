%###################################  计算图像相干性   ##################################%
%###################################  Arthur:Unknown  ##################################%
%################################### date：2017.6.28  ##################################%
%————————————————————————————————————————————————%
%——————————————————————更新状态——————————————————————%
% 版 本 号：1.0.0
% 近期修改：无
% 执 行 人：李莫凡
% 审 核 人：无
%————————————————————————————————————————————————%
%——————————————————————函数介绍——————————————————————%
% 函数功能：
% 计算两幅图像的相干系数，输出复相干系数图
%————————————————————————————————————————————————%
%                    参数名             参数定义              数据尺寸          数据类型
% 输入参数：         master              主图像            二维矩阵（M*N）       double
%                   slave_reg           副图像            二维矩阵（M*N）       double
%             search_block_col\row  计算窗横\纵向尺寸      单个元素（1*1）       double
%————————————————————————————————————————————————%
%                    参数名             参数定义              数据尺寸          数据类型
% 输出参数：     correlation_map       相关系数图          二维矩阵（M*N）       double
%————————————————————————————————————————————————%
function [corr_poly_C]=Correlation_Cal(master,slave_reg)
%%
% ?üD?ê±??: 2017.6.28
% ?????à?é?μêy
% ×￠òa : 
%            ó?óú×÷í? 2?±￡′?êy?Y
%
win_size_corr =4;
[Na,Nr]=size(master);
master_corr = zeros(Na+win_size_corr,Nr+win_size_corr);
slave_corr = zeros(Na+win_size_corr,Nr+win_size_corr);
master_corr(win_size_corr/2+1:end-win_size_corr/2,win_size_corr/2+1:end-win_size_corr/2) = master;
slave_corr(win_size_corr/2+1:end-win_size_corr/2,win_size_corr/2+1:end-win_size_corr/2) = slave_reg;
corr_poly_C = zeros(Na,Nr);
corr_poly_R = zeros(Na,Nr);
for ii = 1 + win_size_corr/2 : Na + win_size_corr/2
%     ii
    for jj = 1 + win_size_corr/2: Nr + win_size_corr/2
        corr_master = master_corr(ii-win_size_corr/2:ii+win_size_corr/2,jj-win_size_corr/2:jj+win_size_corr/2);
        corr_slave = slave_corr(ii-win_size_corr/2:ii+win_size_corr/2,jj-win_size_corr/2:jj+win_size_corr/2);
        %?????′?￠êμ?à?é?μêy
        cc_1 = corr_master.*conj(corr_slave);
        cc_2 = corr_master.*conj(corr_master);
        cc_3 = corr_slave.*conj(corr_slave);
        rr_1 = abs(corr_master).*abs(corr_slave);
        rr_2 = (abs(corr_master)).^2;
        rr_3 = (abs(corr_slave)).^2;
        corr_poly_C(ii-win_size_corr/2,jj-win_size_corr/2) = abs(sum(cc_1(:)))/sqrt(sum(cc_2(:))*sum(cc_3(:)));
        corr_poly_R(ii-win_size_corr/2,jj-win_size_corr/2) = sum(rr_1(:))/(sqrt(sum(rr_2(:)))*sqrt(sum(rr_3(:))));
    end
end
a=1;
end
