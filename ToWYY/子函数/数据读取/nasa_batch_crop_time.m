% 本函数的作用是按照开始结束时刻批量的截取星历数据
function [Shike_Line_Str_crop_Cell,XYZ_crop_Cell] = nasa_batch_crop_time(dir_path,start_time_vector,end_time_vector)
%% 读取文件夹中所有星历文件并合成一组参数
[Shike_Line_Str,XYZ] = nasa_batch_read(dir_path);
Shike_Line_Str = string(Shike_Line_Str);
xingli = [Shike_Line_Str,XYZ];

%% 按照输入的开始结束时间段向量批量的截取星历
% 计算需要进行的截取次数
num = length(start_time_vector);
% 定义截取后星历的存储空间
Shike_Line_Str_crop_Cell = cell(num,1);
XYZ_crop_Cell = cell(num,1);
% 循环遍历所有需要截取的时间段
for i = 1:num
    % 将时间段开始与结束时间转化为一种字符串的形式
    start_time = datestr(start_time_vector(i,1),'yyyy-mm-ddTHH-MM-SS-FFFZ');
    end_time = datestr(end_time_vector(i,1),'yyyy-mm-ddTHH-MM-SS-FFFZ');
    % 获取需要进行截取的范围
    start_time_place = strcmp(start_time,xingli(:,1));
    start_place_num = find(start_time_place == 1);
    end_time_place = strcmp(end_time,xingli(:,1));
    end_place_num = find(end_time_place == 1);
    % 进行星历的截取
    xingli_crop = xingli(start_place_num:end_place_num,:);
    Shike_Line_Str_crop = xingli_crop(:,1);
    XYZ_crop = str2double(xingli_crop(:,2:4));
    % 将信息装入元胞数组
    Shike_Line_Str_crop_Cell{i,1} = Shike_Line_Str_crop;
    XYZ_crop_Cell{i,1} = XYZ_crop;
end