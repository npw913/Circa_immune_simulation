function H = arrowPlot3D(X, Y, Z, varargin)
%ARROWPLOT3D Plot 3D curve with directional arrow heads.
%   ARROWPLOT3D(X, Y, Z) plots a 3D curve with arrows indicating direction.
%
%   Options:
%       'numArrows'     Number of arrows (default: 2)
%       'color'         Color of curve and arrows (default: [0, 0.447, 0.741])
%       'LineWidth'     Curve line width (default: 0.5)
%       'scale'         Arrow size scaling factor (default: 1)
%       'headSize'      Arrow head size (default: 0.1)
%
%   Example:
%       t = linspace(0, 10*pi, 1000);
%       x = t.*cos(t);
%       y = t.*sin(t);
%       z = t;
%       arrowPlot3D(x, y, z, 'numArrows', 5, 'color', 'r', 'headSize', 0.2)

% 默认参数

defaults = struct(...
    'numArrows', 2,...
    'color', [0, 0.447, 0.741],...
    'LineWidth', 0.5,...
    'scale', 1,...
    'headSize', 0.05,...
    'axis', 'auto',...
    'LineVisible',1,...
    'SampleIndices',-1 ...
);


% 解析用户输入参数
params = parseParams(varargin, defaults);

if size(X,1)>1  % 列向量转为行向量
    X = X';
    Y = Y';
    Z = Z';
end
V = [X;Y;Z];

if params.LineVisible==1
% 绘制主曲线
h = plot3(X, Y, Z, 'Color', params.color, 'LineWidth', params.LineWidth);
axis(params.axis);
hold on;
if nargout == 1
    H = h;
end
else
h = plot3(X, Y, Z,'Visible','off');
axis(params.axis);
hold on;
if nargout == 1
    H = h; 
end
end



% 计算曲线总长度
% total_length = sum(sqrt(diff(X).^2 + diff(Y).^2 + diff(Z).^2));
AMAX = [params.axis(2), params.axis(4), params.axis(6)];
AMIN = [params.axis(1), params.axis(3), params.axis(5)];
V_n = diag(1./(AMAX - AMIN)) * (V - diag(AMIN) * ones(3,length(X)));
X_n = V_n(1,:);
Y_n = V_n(2,:);
Z_n = V_n(3,:);



% 等弧长采样箭头位置
if params.SampleIndices == -1
    sample_indices = round(linspace(1, length(X), params.numArrows+2));
    sample_indices = sample_indices(2:end-1); % 排除首尾点
else
    sample_indices = params.SampleIndices;
end

% 生成基础箭头模型（沿x轴方向）
[arrowX, arrowY, arrowZ] = create3DArrow(params.headSize*params.scale);
DT = delaunay(arrowY, arrowZ);
% 沿曲线放置箭头
for idx = sample_indices
    % 获取当前点位置
    pos = [X_n(idx), Y_n(idx), Z_n(idx)];
    
    % 计算切线方向（归一化）
    if idx == length(X)
        tangent = [X_n(idx)-X_n(idx-1), Y_n(idx)-Y_n(idx-1), Z_n(idx)-Z_n(idx-1)];
    else
        tangent = [X_n(idx+1)-X_n(idx), Y_n(idx+1)-Y_n(idx), Z_n(idx+1)-Z_n(idx)];
    end
    tangent = tangent / norm(tangent);
    
    % 计算旋转矩阵
    R = vecRotation([1 0 0], tangent); % 从x轴旋转到切线方向
    
    % 应用旋转和平移
    rotatedArrow = (R*[arrowX; arrowY; arrowZ]);
    pos = diag(pos)* ones(3, length(arrowX(1,:)));
    finalArrow = rotatedArrow + pos;
    finalArrow = diag(AMIN) * ones(3,length(arrowX(1,:))) + diag((AMAX - AMIN))*finalArrow;

    
    % 绘制箭头
%     fill3(finalArrow(1,:), finalArrow(2,:), finalArrow(3,:), ...
%           params.color, 'EdgeColor', 'none');

    trisurf(DT, finalArrow(1,:), finalArrow(2,:), finalArrow(3,:), 'FaceColor', params.color, 'FaceAlpha', 1, 'EdgeColor', 'none');

end

% axis vis3d;
view(3);
end

% 参数解析子函数
function params = parseParams(inputs, defaults)
params = defaults;
for i = 1:2:length(inputs)
    if isfield(defaults, inputs{i})
        params.(inputs{i}) = inputs{i+1};
    end
end
end

% 三维箭头生成子函数
function [X, Y, Z] = create3DArrow(headSize)
% 生成锥形箭头（沿x轴方向）
theta = linspace(0, 2*pi, 64);
r = headSize/2;
z = r/2 * cos(theta);
y = r/2 * sin(theta);
x_base = zeros(size(z)) - r;

% 组合顶点
% X = [x_base, -headSize*1.5*ones(1,8), -headSize*1.5*ones(1,8)]';
% Y = [y, y*0.8, y*0.8]';
% Z = [z, z*0.8, z*0.8]';
X = x_base';
Y = y';
Z = z';
% 添加尖端
X = [X; 0]';
Y = [Y; 0]';
Z = [Z; 0]';

end

% 三维旋转矩阵计算子函数
function R = vecRotation(vec1, vec2)
% 计算从vec1到vec2的旋转矩阵
vec1 = vec1/norm(vec1);
vec2 = vec2/norm(vec2);

axiss = cross(vec1, vec2);
if norm(axiss) < eps
    R = eye(3);
else
    axiss = axiss/norm(axiss);
    angle = acos(dot(vec1, vec2));
    R = axang2rotm([axiss, angle]);
end
end