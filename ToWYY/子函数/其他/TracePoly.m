%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  归一化曲线拟合  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Arthur：LMF   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Date：2021.11.17 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%———————————————————————————————————————%
% 版本：1.0.0
% 最新改动：无
%———————————————————————————————————————%
function [Trace_Res] = TracePoly(Trace_Ori)
% 输入随时间变化的Az或El向量，输出拟合后相同时间点的Az或El向量

    old_x = (1:length(Trace_Ori)).';
    old_y = Trace_Ori;
    % 时间跨度相对于角度变化而言较大
    % 尺度差异较大的横纵坐标会导致拟合失效
    % 需要对样本坐标进行归一化，即(data - min)/(max - min)：
    old_x_max = max(max(old_x)); % 横坐标归一化因子
    old_y_max = max(max(old_y)); % 纵坐标归一化因子
    old_x_min = min(min(old_x)); % 横坐标归一化因子
    old_y_min = min(min(old_y)); % 纵坐标归一化因子
    x_for_interp = (old_x(:) - old_x_min)/(old_x_max - old_x_min);
    y_for_interp = (old_y(:) - old_y_min)/(old_y_max - old_y_min);
    
    p = polyfit(x_for_interp,y_for_interp,4); % 至少需要4阶拟合才能和原轨迹吻合
    new_y = polyval(p,x_for_interp);
    
    Trace_Res = new_y * (old_y_max - old_y_min) + old_y_min; % 重采样纵坐标还原到原始尺度
    
    % 检查效果用
    %{
        figure,plot(Trace_Res,'-*')
        hold on, plot(Trace_Ori,'-o')
    %}
end