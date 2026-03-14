function phase_flat_earth = calculate_Phase_flat(Para)
% 该函数用来计算平地相位
% 输入量：
% 1）R0_RCMC     是随距离线变化的斜距；
% 2）Parameter   代表计算平地相位需要的一些参数，包括：Naz，H，lamda；
% 3）B           代表基线距离；
% 4）theta_B     代表基线倾角；
% 返回值：phase_flat_earth 是所需要的平地相位。
%
% 程序截止到：2014.12.17. 19:57 p.m.

%%
Naz = Para.Naz;
B = Para.B;
lamda = Para.lamda;

for j = 1:length(Para.lat_net)
    for k = 1:length(Para.lon_net)
        r_range_A(j,k) = sqrt((Para.x_radar(round((Naz-B)/2))-Para.flat_x(j,k)).^2 + (Para.y_radar(round((Naz-B)/2))-Para.flat_y(j,k)).^2 + (Para.z_radar(round((Naz-B)/2)) - Para.flat_z(j,k)).^2);
        r_range_B(j,k) = sqrt((Para.x_radar(round((Naz+B)/2))-Para.flat_x(j,k)).^2 + (Para.y_radar(round((Naz+B)/2))-Para.flat_y(j,k)).^2 + (Para.z_radar(round((Naz+B)/2)) - Para.flat_z(j,k)).^2);
    end
end

phase_flat_earth = -4*pi/lamda.*(r_range_A-r_range_B);  % 计算得到的平地相位;
                                                        % 注意，负号很关键！！这由SAR的回波特点决定

end
