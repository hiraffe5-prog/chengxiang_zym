% 本函数的作用是对比配准偏移量真值与该点拟合结果的相对关系
function fit_error_analysis(threshold_cor)
%% 读取拟合前的偏移信息和拟合后的偏移信息
load (['./SaveData/Match_subpix_',num2str(0),'.mat'])           % 读取拟合前的精配准偏移量
load (['./SaveData/Match_subpix_',num2str(0),'_add.mat'])      % 读取拟合后的精配准偏移量                                  

%% 构建相干性掩膜
mask = double(max_cor > threshold_cor);
index = find(mask == 0);
mask(index) = nan;

%% 计算拟合后参考点的偏移量
row_offset_after = fitobjectx(row_pix_master,column_pix_master) - row_pix_master;
col_offset_after = fitobjecty(row_pix_master,column_pix_master) - column_pix_master;

%% 展示拟合前后参考点偏移量差距
figure; imagesc((row_offset_after-row_offset).*mask); title('拟合前后参考点行向偏移量差距')
disp(['行向平均偏移量：',num2str(mean(mean(abs((row_offset_after-row_offset).*mask))))]);
figure; imagesc((col_offset_after-col_offset).*mask); title('拟合前后参考点列向偏移量差距')
disp(['列向平均偏移量：',num2str(mean(mean(abs((col_offset_after-col_offset).*mask))))]);

end