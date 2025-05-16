function plot_sobol_indices(S, ST, param_names)
% 绘制Sobol敏感性指数的分组条形图
% 输入：
%   S          : 一阶敏感性指数，k×1向量
%   ST         : 总效应敏感性指数，k×1向量
%   param_names: 参数名称，cell数组，例如{'X1', 'X2', 'X3'}

    % 确保输入为列向量
    S = S(:);
    ST = ST(:);

    % 按ST从大到小排序，并获取排序索引
    [~, sorted_idx] = sort(ST, 'descend');
    
    % 根据排序索引调整S、ST、param_names
    S = S(sorted_idx);
    ST = ST(sorted_idx);
    param_names = param_names(sorted_idx);
    
    % 合并数据为矩阵，每行对应一个参数，每列为S和ST
    data = [S, ST];
    
    % 创建图形
    figure;
    hBar = bar(data, 'grouped', 'BarWidth', 0.5);
%     hBar.BarWidth = 0.5;
    set(gcf, 'Position', [100 100 800 400]); % 设置图形大小[宽度, 高度]
    
    % 设置颜色（使用MATLAB默认的蓝色和橙色）
    colors = lines(2); % 获取前两种颜色
    hBar(1).FaceColor = colors(1,:); % 第一组条形颜色（S_i）
    hBar(2).FaceColor = colors(2,:); % 第二组条形颜色（ST_i）
    
    % 设置坐标轴标签和标题
    xlabel('Parameters', 'FontSize', 12, 'FontName', 'Arial');
    ylabel('Sensitivity Index', 'FontSize', 12, 'FontName', 'Arial');
    title('Sobol Sensitivity Indices', 'FontSize', 14, 'FontName', 'Arial');
    
    % 设置x轴刻度标签（参数名称）
    xticks(1:length(param_names));
    xticklabels(param_names);
    ax = gca; 
    set(ax, 'FontSize', 12, 'FontName', 'Arial');
    ax.LineWidth = 1; % 加粗边框
    ax.TickLength = [0.008, 0.008]; % 前者设置y轴 tick的长度，后者设置x轴
    set(ax, 'TickDir', 'in')
    
    % 添加图例
    legend({'First-order (S_i)', 'Total-effect (ST_i)'}, ...
           'FontSize', 12, 'Location', 'northeast', 'FontName', 'Arial');
    
    % 设置网格线
    grid on;
    set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3);
    
    % 调整Y轴范围（根据数据自动扩展5%空间）
    ylim([0, 1.2]);
    box on
    
    % 保存为高分辨率图片（PNG/PDF）
%     exportgraphics(gcf, 'Sobol_Sensitivity.png', 'Resolution', 300); % 或替换为'PDF'
end