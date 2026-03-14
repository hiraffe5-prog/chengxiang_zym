function [Shike_Line_Str_crop,XYZ_crop]=nasa_crop_time(Ephemeris_Path,start_time,end_time)
% Shike_Line_Str_crop:N*-,string
% XYZ_crop:N*3，double

% Ephemeris_Path = 'earth_trace_visible.csv';
% start_time = '2021-01-05T09-01-00-000Z';
% end_time = '2021-01-06T17-47-00-000Z';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A = readmatrix(Ephemeris_Path,'OutputType','string');
% xingli = A;
[Shike_Line_Str,XYZ]=nasa_csv_read(Ephemeris_Path);
Shike_Line_Str = string(Shike_Line_Str);
xingli = [Shike_Line_Str,XYZ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 字符串对比
start_time_place = strcmp(start_time,xingli(:,1));
start_place_num = find(start_time_place == 1);
end_time_place = strcmp(end_time,xingli(:,1));
end_place_num = find(end_time_place == 1);
% 星历截取
xingli_crop = xingli(start_place_num:end_place_num,:);
Shike_Line_Str_crop=xingli_crop(:,1);
XYZ_crop=str2double(xingli_crop(:,2:4));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 保存截取后的星历
% delete('xx_xxx_nasa_trace_crop.csv');
% fid = fopen('xx_xxx_nasa_trace_crop.csv','a');
% if (fid < 0)
%     errordlg('File creation failed!','Error');
% end
% 
% fprintf(fid, 'UTC calendar date,X (km),Y (km),Z (km)\n');
% 
% for ii = 1 : size(xingli_crop,1)
%     fprintf(fid, '%s,%s,%s,%s\n',xingli_crop(ii,1),xingli_crop(ii,2),xingli_crop(ii,3),xingli_crop(ii,4));
% end
% 
% fclose(fid);
end