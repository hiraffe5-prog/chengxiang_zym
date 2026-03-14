% 本函数的作用是批量读入文件夹中的星历文件并将所有信息按时间顺序整合到一个文件中，因此文件夹中所有文件覆盖的时间范围之和必须是连续的
function [Shike_Line_Str,XYZ] = nasa_batch_read(dir_path)
%% 获取存储星历的文件夹中所有星历文件的路径
files = dir([dir_path,'*.csv']);            % 所有星历文件的路径
file_num = length(files);                   % 星历文件数

%% 循环读取所有文件中的数据
% 设置时间与星历信息的存储元胞数组
Shike_Line_Str_Cell = cell(file_num,1);
XYZ_Cell = cell(file_num,1);
% 设置记录各个文件起始时间点与终止时间点的向量
file_start_time = zeros(file_num,1);
file_end_time = zeros(file_num,1);
% 循环读取所有数据
for i = 1:file_num
    [Shike_Line_Str_tmp,XYZ_tmp] = nasa_csv_read(fullfile(files(i).folder,files(i).name));                                   % 读入第i个星历文件的数据
    file_start_time(i,1) = datenum(Shike_Line_Str_tmp(1,:),'yyyy-mm-ddTHH-MM-SS-FFFZ');         % 保存读入数据的起始时刻（日期序列的形式）
    file_end_time(i,1) = datenum(Shike_Line_Str_tmp(end,:),'yyyy-mm-ddTHH-MM-SS-FFFZ');         % 保存读入数据的结束时刻（日期序列的形式）
    Shike_Line_Str_Cell{i,1} = Shike_Line_Str_tmp;                                              % 将读入的数据存入元胞数组中
    XYZ_Cell{i,1} = XYZ_tmp;
end

%% 对文件中的数据进行拼接
% 对文件数据的开始时间进行排序
[~,flag] = sort(file_start_time);
% 按照排序结果进行数据拼接
for i = 1:length(flag)
    flag_tmp = flag(i);                                                                         % 当前轮到的文件标号
    if i == 1
        Shike_Line_Str = Shike_Line_Str_Cell{flag_tmp,1};
        XYZ = XYZ_Cell{flag_tmp,1};
    else
        if round((file_start_time(flag(i),1) - file_end_time(flag(i-1),1))*24*60) ~= 1          % 判断读取是否符合预期
            error('文件夹中的星历文件不连续!');
        end
        Shike_Line_Str = [Shike_Line_Str;Shike_Line_Str_Cell{flag_tmp,1}];
        XYZ = [XYZ;XYZ_Cell{flag_tmp,1}];
    end
end