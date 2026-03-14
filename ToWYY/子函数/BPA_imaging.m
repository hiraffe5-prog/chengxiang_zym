function [s_ac] = BPA_imaging(Para,raw_data_type)
%% 载入参数
fs = Para.fs;
Naz = Para.Naz;
Nr = Para.Nr;
Nrg = Para.Nrg;
B = Para.B;
Kr = Para.Kr;
Tr = Para.Tr;
R0 = Para.R0;
c = Para.c;
lamda = Para.lamda;

if raw_data_type == 1
    s_echo = Para.s_echoA(1:Naz-B,:);           % 载入天线A对应的原始数据
    x_radar = Para.x_radar(1:Naz-B);
    y_radar = Para.y_radar(1:Naz-B);
    z_radar = Para.z_radar(1:Naz-B);
end
if raw_data_type == 2
    s_echo = Para.s_echoA(B+1:end,:);           % 载入天线B对应的原始数据
    x_radar = Para.x_radar(B+1:end);
    y_radar = Para.y_radar(B+1:end);
    z_radar = Para.z_radar(B+1:end);
end

%% 脉压
% S_range = fft(s_echo,[],2);     % 进行距离向傅里叶变换，零频在两端。
% t_ref = ( -Nr/2 : (Nr/2-1) )/fs;    % 用来生成距离MF的距离时间轴，长度为Nr（Nr<Nrg）；
% t_ref_mtx = ones(Naz-B,1)*t_ref;
% w_range_ref = (abs(t_ref_mtx)) <= ((Tr/2).*ones(Naz-B,Nr));
% s_ref =w_range_ref.* exp(1j*pi*Kr.*(t_ref_mtx).^2); % 复制（发射）脉冲，未加窗。
% 
% 
% s_ref = [s_ref,zeros(Naz-B,Nrg-Nr)];      % 对复制脉冲，末端补零至长度 Nrg；
% s_ref = circshift(s_ref,[0 -Nr/2]);     % 向左循环移位（Nr/2）；
% 
% 
% S_ref = fft(s_ref,[],2);            % 复制脉冲的距离傅里叶变换，零频在两端。
% H_range = conj(S_ref);                % 距离向匹配滤波器，零频在两端。
% % 对距离频谱进行：距离压缩 + 包络校正 + 一次运动补偿
% S_range_c = S_range.*H_range;
% s_rc = ifft(S_range_c,Nrg,2);            % 完成距离压缩，包络校正和一次运补，回到二维时域。
% S_rc = fft(s_rc);
% S_rc = [S_rc(Naz-B,1:ceil(Nrg/2)),zeros(Naz-B,7*Nrg),S_rc(Naz-B,ceil(Nrg/2)+1:end)];
% s_rc =ifft(S_rc);

ref_t=(-Nr/2:Para.Nr/2-1)/Para.Nrg*Para.PRT;       % -Tp/2:Tp/2的时间
ref=exp(1i*pi*Para.Kr*ref_t.^2);                   % 参考信号s：exp(j pi Kr t^2)
ref=conj(fft(ref,Para.Nrg,2));                     % 参考信号频域并取共轭
%fr=[0:Nr/2-1 -Nr/2:-1]/Nr*fs;
%ref=exp(1i*pi*fr.^2/k);
ym=bsxfun(@times,fft(s_echo,Para.Nrg,2),ref);      % 脉压
s_rc=ifft(ym,Para.Nrg,2);                          % 逆fft回到时域
s_rc = circshift(s_rc, [0 Para.Nr/2]);             % 脉压后的信号循环右移，-Tr/2~Tr/2特有，目的是补掉脉压导致的Tr/2时延落后
% S_rc = fft(s_rc,[],2);
% S_rc = [S_rc(:,1:ceil(Nrg/2)),zeros(Naz-B,1*Nrg),S_rc(:,ceil(Nrg/2)+1:end)];
% s_rc =ifft(S_rc,[],2);
%% BP
tr_RCMC = Para.tr+Para.amb_num*Para.PRT;            % 含模糊数的快时间轴，由于计算回波所处采样点数
f_back = zeros(length(Para.lon_net),length(Para.lat_net));             %存储反投结果
for i = 1:Naz-B     % 遍历方位点数
    s_temp = s_rc(i,:);     % 当前方位向上的回波
    rp = sqrt((x_radar(i)-Para.flat_x).^2 + ...
        (y_radar(i)-Para.flat_y).^2 + ...
        (z_radar(i)-Para.flat_z).^2);    %当前方位向雷达与场景各点斜距计算
    n_round=round((rp*2/c*Para.Nrg/Para.PRT-tr_RCMC(1)*Para.Nrg/Para.PRT)+1);   % 当前方位向雷达到场景各点的采样点数计算
    n_round(n_round<1) = 1; n_round(n_round>Nrg) = Nrg;     %避免超范围
    f_back = f_back + s_temp(n_round).*exp(1j*4*pi*rp/lamda);   % 相位补偿，进行各方位向信号的同相叠加
    disp(i/(Naz-B));    
end
% 保相
R_jk = sqrt((x_radar(round((Naz-B)/2))-Para.flat_x).^2 + (y_radar(round((Naz-B)/2))-Para.flat_y).^2 + (z_radar(round((Naz-B)/2)) - Para.flat_z).^2);
f_back_new = f_back.*exp(-1j*4*pi*R_jk/lamda);

s_ac = f_back_new;
% 作图
% 图7——成像结果
figure;
imagesc(Para.lon_net,Para.lat_net,abs(f_back_new));

xlabel('经度');
ylabel('纬度');
if raw_data_type == 1
    title('天线 A，成像结果');
    save('imag_A_points','s_ac');
end
if raw_data_type == 2       
    title('天线 B，成像结果');
    save('imag_B_points','s_ac');
end


end





