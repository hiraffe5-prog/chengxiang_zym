function [X,Y,Z,X0,Y0,Z0] = AER2MCMF(Az,El,Range,B,L,H)
%######################## 将AER坐标转换为MCMF坐标 ########################%
   %% AER转站心坐标系
    x = Range .* cosd(El) .* sind(Az);     % 正东方向
    y = Range .* cosd(El) .* cosd(Az);     % 正北方向
    z = Range .* sind(El);                 % 天顶方向
    %% 站心坐标系转换为月固坐标系
    %
    % 转换矩阵
    X1 = (-sind(L).*x) +  (-sind(B).*cosd(L).*y)  + (cosd(B).*cosd(L).*z);
    Y1 = ( cosd(L).*x)  + (-sind(B).*sind(L).*y)  + (cosd(B).*sind(L).*z);
    Z1 =       0       +  (cosd(B).*y)            + (sind(B).*z);
    % 偏移矩阵
    a = 1738.1e3;                      %  月球椭圆的长半轴
    b = 1736e3;                        %  月球椭圆的短半轴
    e = sqrt(a.^2 - b.^2)/a;           %  月球椭圆的第一偏心率;
    N = a/(1 - e^2*(sind(B)^2));
    X0 = (N + H)*cosd(L)*cosd(B);
    Y0 = (N + H)*sind(L)*cosd(B);
    Z0 = (N*(1 - e.^2) + H)*sind(B);
    % 月固坐标
    X = X0 + X1;
    Y = Y0 + Y1;
    Z = Z0 + Z1;
end