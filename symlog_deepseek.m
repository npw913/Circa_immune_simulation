function symlog_deepseek(varargin)
% SYMLOG 双对称对数坐标轴缩放（支持自定义显示范围）
%   SYMLOG(ax, var, C, 'XLim', [xmin,xmax], ...) 
%   参数说明：
%     - ax: 坐标轴句柄（可选，默认gca）
%     - var: 坐标轴方向（'x','y','z'，默认'xyz'）
%     - C: 缩放常数（默认0）
%     - 'XLim','YLim','ZLim': 指定原始数据范围（例如'YLim',[0,1000]）

    % 默认参数
    ax = [];
    var = 'xyz';
    C = 0;
    lim = struct('x',[], 'y',[], 'z',[]);
    
    % 解析输入参数
    args = varargin;
    while ~isempty(args)
        arg = args{1};
        if isempty(arg)
            args(1) = [];
            continue;
        end
        switch class(arg)
            case 'matlab.graphics.axis.Axes'
                ax = arg;
                args(1) = [];
            case 'char'
                if any(strcmpi(arg, {'xlim','ylim','zlim'}))
                    % 处理范围参数（例如'YLim', [0,1000]）
                    axis_var = lower(arg(1)); % 提取x/y/z
                    lim.(axis_var) = args{2};
                    args(1:2) = [];
                else
                    var = lower(arg);
                    args(1) = [];
                end
            case {'double','single'}
                C = arg;
                args(1) = [];
            otherwise
                error('无效输入类型: %s', class(arg));
        end
    end
    
    % 设置默认坐标轴
    if isempty(ax), ax = gca; end
    
    % 递归处理多轴
    if length(var) > 1
        for v = var
            symlog_deepseek(ax, v, C, lim);
        end
        return
    end
    
    % 核心逻辑
    C = 10.^C; % 转换为实际缩放因子
    
    % 1. 检查历史变换并更新UserData
    userdata = get(ax, 'UserData');
    if isfield(userdata, 'symlog') && isfield(userdata.symlog, var)
        lastC = userdata.symlog.(var);
    else
        lastC = [];
    end
    userdata.symlog.(var) = C;
    set(ax, 'UserData', userdata);
    
    % 2. 转换图形对象并获取极值
    [minVal, maxVal] = transform_objects(ax, var, C, lastC);
    
    % 3. 处理用户指定的范围
    user_lim = lim.(var);
    if ~isempty(user_lim)
        % 将用户指定的原始范围转换为变换后的坐标
        lim_transformed = sign(user_lim) .* log10(1 + abs(user_lim)/C);
        set(ax, [var, 'Lim'], lim_transformed);
    elseif strcmpi(get(ax, [var, 'LimMode']), 'auto')
        set(ax, [var, 'Lim'], [minVal, maxVal]);
    end
    
    % 4. 生成刻度标签
    generate_ticks(ax, var, C, lim.(var));
end

%% 子函数1: 转换图形对象并返回极值
function [minVal, maxVal] = transform_objects(ax, var, C, lastC)
    minVal = inf;
    maxVal = -inf;
    
    % 处理所有图形对象类型
    types = {'line', 'patch', 'rectangle'};
    for t = 1:length(types)
        objs = findobj(ax, 'Type', types{t});
        for ii = 1:length(objs)
            switch types{t}
                case {'line', 'patch'}
                    % 处理线条/面片数据
                    x = get(objs(ii), [var, 'Data']);
                    if ~isempty(lastC)
                        x = sign(x).*lastC.*(10.^abs(x)-1); % 撤销旧变换
                    end
                    x_trans = sign(x).*log10(1 + abs(x)/C); % 应用新变换
                    set(objs(ii), [var, 'Data'], x_trans);
                    
                case 'rectangle'
                    % 处理矩形位置
                    pos = get(objs(ii), 'Position');
                    switch var
                        case 'x'
                            x = [pos(1), pos(1)+pos(3)];
                        case 'y'
                            x = [pos(2), pos(2)+pos(4)];
                    end
                    if ~isempty(lastC)
                        x = sign(x).*lastC.*(10.^abs(x)-1);
                    end
                    x_trans = sign(x).*log10(1 + abs(x)/C);
                    
                    % 更新位置
                    switch var
                        case 'x'
                            pos(1) = x_trans(1);
                            pos(3) = x_trans(2) - x_trans(1);
                        case 'y'
                            pos(2) = x_trans(1);
                            pos(4) = x_trans(2) - x_trans(1);
                    end
                    set(objs(ii), 'Position', pos);
                    x_trans = x_trans(:); % 统一为列向量
            end
            
            % 更新极值
            if ~isempty(x_trans)
                minVal = min(minVal, min(x_trans));
                maxVal = max(maxVal, max(x_trans));
            end
        end
    end
    
    % 处理无数据情况
    if isinf(minVal), minVal = 0; maxVal = 0; end
end

%% 子函数2: 生成刻度标签（仅显示在用户指定范围内）
function generate_ticks(ax, var, C, user_lim)
    % 原始刻度生成逻辑
    t0 = ceil(log10(C)) : ceil(log10(abs(user_lim(2))/C*10));
    t1 = 10.^t0;
    
    % 合并正负刻度
    t_all = [ -fliplr(t1), 0, t1 ];
    lbl_all = arrayfun(@(x)sprintf('10^{%d}',x), [ -fliplr(t0), 0, t0 ], 'UniformOutput',false);
    lbl_all{t_all==0} = '0';
    
    % 转换到symlog空间
    t_trans = sign(t_all).*log10(1 + abs(t_all)/C);
    
    % 过滤用户范围
    if ~isempty(user_lim)
        lim_trans = sign(user_lim).*log10(1 + abs(user_lim)/C);
        valid = (t_trans >= lim_trans(1)) & (t_trans <= lim_trans(2));
        t_trans = t_trans(valid);
        lbl_all = lbl_all(valid);
    end
    
    % 设置主刻度和标签
    set(ax, [var, 'Tick'], t_trans, [var, 'TickLabel'], lbl_all);
    
    % 生成次刻度
    s_all = [];
    for t = t_all
        if t == 0
            continue;
        elseif t > 0
            s = (2:9) * t;
        else
            s = (2:9) * abs(t) * (-1);
        end
        s_all = [s_all, s];
    end
    
    % 转换到symlog空间
    s_trans = sign(s_all) .* log10(1 + abs(s_all)/C);
    
    % 过滤到用户范围
    if ~isempty(user_lim)
        valid_s = (s_trans >= lim_trans(1)) & (s_trans <= lim_trans(2));
        s_trans = s_trans(valid_s);
    end
    
    % 设置次刻度（兼容MATLAB版本）
    var_axis = [upper(var), 'Axis'];  % 关键修复：确保属性名大写（如'YAxis'）
    if isprop(ax, var_axis)
        ax.(var_axis).MinorTickValues = s_trans;
    else
        warning('当前MATLAB版本不支持直接设置次刻度位置。请升级至R2015b或更高版本。');
    end
    
    % 启用次刻度显示
    set(ax, [var, 'MinorTick'], 'on');
end

