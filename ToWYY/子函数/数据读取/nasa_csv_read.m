function [Shike_Line_Str,XYZ] = nasa_csv_read(EARTH_emp_path)
% Shike_Line_Str:N*-,string
% XYZ:N*3,string

%EARTH_emp_path = 'IAU_MARS_DE421_2021_12_10T04_00_00_000_2021_12_12T03_59_00_000.csv';
trace_E = readmatrix(EARTH_emp_path,'OutputType','string');
trace_E = trace_E(13:end-4,1:3);
Num_E = size(trace_E,1);
earth_trace = zeros(Num_E,4);
earth_trace = string(earth_trace);
for ii = 1:Num_E
    t_3 = trace_E(ii,1);
    t_3 = regexp(t_3,',', 'split');
    earth_trace(ii,2) = str2double(t_3(4));
    earth_trace(ii,3) = str2double(t_3(5));
    earth_trace(ii,4) = str2double(t_3(6));
    t_e_3 = regexp(t_3(1),'[TZ:. ]', 'split');
    t_mm = str2double(t_e_3(5));
    t_mm = sprintf('%03d',floor(t_mm/100));
    earth_trace(ii,1) = [char(t_e_3(1)),'T',char(t_e_3(2)),'-',char(t_e_3(3)),'-',char(t_e_3(4)),'-',t_mm,'Z'];
end

Shike_Line_Str = earth_trace(:,1);
XYZ = str2double(earth_trace(:,[2,3,4]));
end