function demo_fig8
Time_circadian_clock = xlsread('Timing.xlsx'); 
Par_new = [145 / 8 ...% gamma_PTLN
    1.573036637...% c_PTLN_cell
    246.434246236928 ...% K_CCL19_DC
    505.3267033...% K_CCL21_DC
    0.502396976147 ...% alpha
    1.332474469 ...% d_DC
    1000 ...% a_0_CCL
    5389.733217150947...% a_DC_CCL
    3.534061671613...% K_DC_CCL
    7.408]; % d_CCL
colors = lines(100);


y0 = [47.5877748302043 0 0 0 0 0 0]; 
T70 = 0 : 0.001 : 7 / 24; 
T190 = 0 : 0.001 : 19 / 24; 
T = 0 : 0.001 : 4; 
[~, Y70] = ode15s(@(t, y)ode_DCCCL_new_ut(t, y, Par_new, Time_circadian_clock, 0), T70, y0); 
[~, Y190] = ode15s(@(t, y)ode_DCCCL_new_ut(t, y, Par_new, Time_circadian_clock, 0), T190, y0); 
y70 = Y70(end, :); 
y190 = Y190(end, :); 

% BMDC
y70(4) = 10 ^ 6 / 145000; 
y190(4) = 10 ^ 6 / 145000; 

[~, Y7] = ode15s(@(t, y)ode_DCCCL_new_ut(t, y, Par_new, Time_circadian_clock, 7), T, y70); 
[~, Y19] = ode15s(@(t, y)ode_DCCCL_new_ut(t, y, Par_new, Time_circadian_clock, 19), T, y190); 
ut_ZT7 = zeros(1, length(T)); 
ut_ZT19 = zeros(1, length(T)); 
for i = 1 : length(T)
    dydt = ode_DCCCL_new_ut(T(i), Y7(i, :), Par_new, Time_circadian_clock, 7); 
    ut_ZT7(i) = dydt(5); 
    dydt = ode_DCCCL_new_ut(T(i), Y19(i, :), Par_new, Time_circadian_clock, 19); 
    ut_ZT19(i) = dydt(5); 
end

UT = 0 : 0.01 : 7; 
DCPT = 0 : 0.01 : 7; 
[X, Y] = meshgrid(DCPT, UT); 
% Z = zeros(length(UT),length(DCPT),3);  % CCL19/21_DCs
% for i = 1:length(DCPT)
%     for j = 1:length(UT)
%         DC = 0:0.001:150;
%         bx =Par_new(1)*Par_new(2)*Par_new(5)*DCPT(i)/Par_new(6);
%         bxut =UT(j)/Par_new(6);
%         k=1;
%             for jj = 1:length(DC)-1
%                  CCL_jj =  by*DC(jj).^2./(Kd^2+DC(jj).^2);
%                  CCL_jplus1 = by*DC(jj+1).^2./(Kd^2+DC(jj+1).^2);
%                  a = DC(jj) -  bx * CCL_jj.^2./(Kc^2+CCL_jj.^2) - bxut ;
%                  b = DC(jj+1) -  bx * CCL_jplus1.^2./(Kc^2+CCL_jplus1.^2) - bxut;
%                  if (a*b<0)
%                      Z(j,i,k)=(CCL_jj+CCL_jplus1)/2;
%                      k=k+1;
%                  elseif a==0
%                      Z(j,i,k)=CCL_jj;
%                      k=k+1;
%                  else
%                  end
%
%             end
%     end
% end
% Z1 = Z;
% save two_dimension_bifurcation.mat Z1

% Fix points solution(sCCL21_steady_state=sCCL21_ss(DC_PT,ut))
load Z.mat Z
% % 缩小
% % shr = 1:10:701;
% % Z = Z(shr,shr,:);
% 
% % 提取有效离散点
[UTindex_bi,DCPTindex_bi] = find(Z(:,:,2));%在saddle node layer寻找非零元素，确定bistable的区域
index2 = find(Z(:, :, 2)); %在saddle point layer寻找非零元素，确定bistable的区域
index1 = find(Z(:, :, 1)); 
index3 = find(Z(:, :, 3)); 
% 
[row,col]=find(Z(:, :, 2));  % 行代表ut（y轴），列代表DCpt（x轴）
% 
% 
% 
l_Z = Z(:, :, 1); 
s_Z = Z(:, :, 2); 
h_Z = Z(:, :, 3); 
spl = l_Z(index1); 
% spl_true = l_Z(index2); 
sps = s_Z(index2); 
sph = h_Z(index3); 
% a=s_Z(index);
% a=Z(UTindex_bi,DCPTindex_bi,2);
f=figure("Name",'Fig8',NumberTitle='off');
f.Color = [1 1 1];
% s=surf(X(UTindex_bi,DCPTindex_bi),Y(UTindex_bi,DCPTindex_bi),Z(UTindex_bi,DCPTindex_bi,2),'FaceAlpha',0.1);
% s.EdgeColor = 'none';
% s=surf(X,Y,Z(:,:,1),'FaceAlpha',0);
% s.EdgeColor='none';
% hold on
% s1=surf(X,Y,Z(:,:,2),'FaceAlpha',0.1);
% s1.EdgeColor='none';
% s2=surf(X,Y,Z(:,:,3),'FaceAlpha',0.1);
% s2.EdgeColor='none';
C = 2; 
TT = tiledlayout(2, 2); 
% spl = symlog10(spl,C);
% sps = symlog10(sps,C);
% sph = symlog10(sph,C);
% Z_ZT7 = symlog10(Y7(:,2),C);
% Z_ZT19 = symlog10(Y19(:,2),C);
j = nexttile(2); 
XX_s = [X(index1); X(index2); X(index3)]; 
YY_s = [Y(index1); Y(index2); Y(index3)]; 
ZZ_s = [spl; sps; sph]; 
% griddata =
%     p1 = plot3(Y7(:, 4), ut_ZT5, Y7(:, 2), LineWidth = 1.5, Color = '#1F77B4'); 
grid on
hold on
axis square
% p2 = plot3(Y19(:, 4), ut_ZT17, Y19(:, 2), LineWidth = 1.5, Color = '#D62728'); 
% p3 = scatter3(X(index1), Y(index1), spl, 0.5, '.', 'MarkerFaceAlpha', 0.1, 'MarkerFaceColor', colors(1, :), 'MarkerEdgeColor', colors(1, :), 'MarkerEdgeAlpha', 0.1); 

% p4 = scatter3(X(index2), Y(index2), sps, 0.5, '.', 'MarkerFaceAlpha', 0.1, 'MarkerFaceColor', colors(2, :), 'MarkerEdgeColor', colors(2, :), 'MarkerEdgeAlpha', 0.1); 
% scatter3(X(index3), Y(index3), sph, 0.5, '.', 'MarkerFaceAlpha', 0.1, 'MarkerFaceColor', colors(1, :), 'MarkerEdgeColor', colors(1, :), 'MarkerEdgeAlpha', 0.1);
% scatter(X(index3), Y(index3),'ro');
shp1 = alphaShape([0 7 7 0]',[0 0 7 7]',100);
p1=plot(shp1,EdgeColor="none",FaceColor=colors(1,:),FaceAlpha=0.3);
hold on
shp = alphaShape(X(index3),Y(index3),0.1);
p2=plot(shp,EdgeColor="none",FaceColor=colors(2,:),FaceAlpha=0.3);
p3=arrowPlot(Y7(:, 4), ut_ZT7,...
        'number', 1, 'color', [0.122 0.467 0.706], 'LineWidth', 1.5, 'scale', 1,'limit',[0 7 0 7],'size',8);
p4=arrowPlot(Y19(:, 4), ut_ZT19,...
        'number', 1, 'color', [0.839 0.153 0.157], 'LineWidth', 1.5, 'scale', 1,'limit',[0 7 0 7],'size',8);

xlabel('DC_{PT} (10^3/\mul)', FontSize = 13); 
% ylabel('u(10^3\cdot\mul^{-1}\cdotd^{-1})',FontSize=13);
ylabel('u (10^3/(\mul·day))', FontSize = 13); 
% zlabel('sCCL19/sCCL21(pg/ml)', FontSize = 13); 
    set(gca, 'FontName', 'Arial'); 
    ax = gca; 
    ax.LineWidth = 1; % 加粗边框
    ax.TickLength = [0.008, 0.008]; % 前者设置y轴 tick的长度，后者设置x轴
    axis([0 7 0 5.5]);
    axis square
    
    legend([p1, p2, p3, p4], {'Monostable', 'Bistable','ZT7','ZT19'},'Location','Northwest'); 
% plot(DCPT(DCPTindex_bi),UT(UTindex_bi),'r.');
% hold off
nexttile(1)

%三角剖分前做区域切分
xmin_s = min(X(index2)); 
xmax_s = max(X(index2)); 
ymax_s = max(Y(index2)); 
s_X = X(index2); 
s_Y = Y(index2); 
l_X = X(index1); 
l_Y = Y(index1); 
h_X = X(index3); 
h_Y = Y(index3); 

s_X = round(s_X*100)*1e-02;
s_Y = round(s_Y*100)*1e-02;
l_X = round(l_X*100)*1e-02;
l_Y = round(l_Y*100)*1e-02;
h_X = round(h_X*100)*1e-02;
h_Y = round(h_Y*100)*1e-02;
N = 1; 
for i = 1:N
    %         p_command_string_5 = "p"+num2str(i)+"_"+num2str(j)+".EdgeColor = 'none';";
%         eval(p_command_string_1);

    com_str = "index_s"+num2str(i)+"= find((s_X>(xmin_s+(i-1)*(xmax_s-xmin_s)/N))&(s_X<(xmin_s+i*(xmax_s-xmin_s)/N)));";
    eval(com_str);
    com_str = "DT_s"+num2str(i)+" = delaunay(s_X(index_s"+num2str(i)+"),s_Y(index_s"+num2str(i)+"));";
    eval(com_str);
%     com_str = "trisurf(DT_s"+num2str(i)+",s_X(index_s"+num2str(i)+"),s_Y(index_s"+num2str(i)+"),sps(index_s"+num2str(i)+...
%         "),'FaceColor',colors(2,:),'FaceAlpha',0.2,'EdgeColor',colors(2,:));";
%     eval(com_str);
%     hold on
%     com_str = "trisurf(DT_s"+num2str(i)+",X(index_s"+num2str(i)+"),Y(index_s"+num2str(i)+"),l_Z(index_s"+num2str(i)+...
%         "),'FaceColor',colors(1,:),'FaceAlpha',0.2,'EdgeColor','none');";
%     eval(com_str);
end
DT_sl = delaunay(l_X, l_Y); 
DT_ss = delaunay(s_X, s_Y); 
% DT_sh = delaunay(h_X,h_Y);
% err_i = find((s_X==0.28)&(s_Y==1.32));
% index_DT_s1 = find(DT_ss==err_i);
% DT_ss(index_DT_s1)=0;

% 削
index_row = []; 
for i = 1 : length(DT_ss(:, 1))
    i1 = DT_ss(i, 1); 
    i2 = DT_ss(i, 2); 
    i3 = DT_ss(i, 3); 
    
    sum_abs_z = abs(sps(i1) - sps(i2)) + abs(sps(i2) - sps(i3)) + abs(sps(i3) - sps(i1)); 
%     if (DT_ss(i,1)*DT_ss(i,2)*DT_ss(i,3))==0
    if sum_abs_z > 95
        index_row = [index_row; i]; 
    end
end
DT_ss(index_row, :) = []; 
index_row = []; 
for i = 1 : length(DT_sl(:, 1))
    i1 = DT_sl(i, 1); 
    i2 = DT_sl(i, 2); 
    i3 = DT_sl(i, 3); 
    sum_abs_z = abs(spl(i1) - spl(i2)) + abs(spl(i2) - spl(i3)) + abs(spl(i3) - spl(i1)); 
    if sum_abs_z > 350
        index_row = [index_row; i]; 
    end
end
DT_sl(index_row, :) = []; 

% 补
l_h_X_h = 0.28 : 0.01 : 7; 
l_h_X_h = l_h_X_h'; 
l_h_Y_h = zeros(length(l_h_X_h), 1); 
l_h_Z_h = zeros(length(l_h_X_h), 1); 
for i = 29 : 701
index=find(col==i);
index=max(row(index));
l_h_Y_h(i-28)=0.01*(index-1);
l_h_Z_h(i-28)=h_Z(round(100*(l_h_Y_h(i-28)+0.01)),round(100*(l_h_X_h(i-28)+0.01)));
end

l_h_X_l = 0.28: 0.01 : 7; 
l_h_X_l = l_h_X_l'; 
l_h_Y_l = l_h_Y_h+0.01; 
l_h_Z_l = zeros(length(l_h_X_l), 1); 

s_l_X_l = 0.28: 0.01 : 7; 
s_l_X_l = s_l_X_l'; 
s_l_Y_l = l_h_Y_h;
s_l_Z_l = zeros(length(s_l_X_l), 1);

s_l_X_s = 0.28: 0.01 : 7; 
s_l_X_s = s_l_X_s';
s_l_Y_s = l_h_Y_h;
s_l_Z_s = zeros(length(s_l_X_l), 1);
for i = 29:701
    l_h_Z_l(i-28)=l_Z(round(100*(l_h_Y_l(i-28)+0.01)),round(100*(l_h_X_l(i-28)+0.01)));  
    s_l_Z_l(i-28)=l_Z(round(100*(s_l_Y_l(i-28)+0.01)),round(100*(s_l_X_l(i-28)+0.01))); 
    s_l_Z_s(i-28)=s_Z(round(100*(s_l_Y_s(i-28)+0.01)),round(100*(s_l_X_s(i-28)+0.01)));
end

h_s_Y_s = ([1:130,132,133]-1)*0.01; 
h_s_Y_s = h_s_Y_s';
h_s_X_s = zeros(length(h_s_Y_s), 1);
h_s_Z_s = zeros(length(h_s_Y_s), 1);

h_s_Y_h = ([1:130,132,133]-1)*0.01; 
h_s_Y_h = h_s_Y_h';
h_s_X_h = zeros(length(h_s_Y_h), 1);
h_s_Z_h = zeros(length(h_s_Y_h), 1);
for j = [1:130,132,133]
%     index=find(col==i);
% index=max(row(index));
% l_h_Y_h(i-28)=0.01*(index-1);
% l_h_Z_h(i-28)=h_Z(round(100*(l_h_Y_h(i-28)+0.01)),round(100*(l_h_X_h(i-28)+0.01)));
    index = find(row==j);
    index = min(col(index));
    if j<131
    h_s_X_s(j)=0.01*(index-1);
    h_s_X_h(j)=0.01*(index-1);
    h_s_Z_s(j)=s_Z(round(100*(h_s_Y_s(j)+0.01)),round(100*(h_s_X_s(j)+0.01)));
    h_s_Z_h(j)=h_Z(round(100*(h_s_Y_h(j)+0.01)),round(100*(h_s_X_h(j)+0.01)));
    else 
    h_s_X_s(j-1)=0.01*(index-1);
    h_s_X_h(j-1)=0.01*(index-1);
    h_s_Z_s(j-1)=s_Z(round(100*(h_s_Y_s(j-1)+0.01)),round(100*(h_s_X_s(j-1)+0.01)));
    h_s_Z_h(j-1)=h_Z(round(100*(h_s_Y_h(j-1)+0.01)),round(100*(h_s_X_h(j-1)+0.01)));        
    end
end

l_h_X = [l_h_X_l;l_h_X_h];
l_h_Y = [l_h_Y_l;l_h_Y_h];
l_h_Z = [l_h_Z_l;l_h_Z_h];

s_l_X = [s_l_X_s;s_l_X_l];
s_l_Y = [s_l_Y_s;s_l_Y_l];
s_l_Z = [s_l_Z_s;s_l_Z_l];

h_s_X = [h_s_X_h;h_s_X_s];
h_s_Y = [h_s_Y_h;h_s_Y_s];
h_s_Z = [h_s_Z_h;h_s_Z_s];

DT_l_h = delaunay(l_h_X,l_h_Y);
DT_s_l = delaunay(s_l_X,s_l_Z);
DT_h_s = delaunay(h_s_Y,h_s_Z);

index_row = []; 
for i = 1 : length(DT_l_h(:, 1))
    i1 = DT_l_h(i, 1); 
    i2 = DT_l_h(i, 2); 
    i3 = DT_l_h(i, 3); 
    sum_abs_z = abs(l_h_Z(i1) - l_h_Z(i2)) + abs(l_h_Z(i2) - l_h_Z(i3)) + abs(l_h_Z(i3) - l_h_Z(i1)); 
    if sum_abs_z > 400
        index_row = [index_row; i]; 
    end
end
DT_l_h(index_row, :) = []; 



index_row = [];
for i = 1 : length(DT_s_l(:, 1))
    i1 = DT_s_l(i, 1); 
    i2 = DT_s_l(i, 2); 
    i3 = DT_s_l(i, 3); 
    sum_abs_x = abs(s_l_X(i1) - s_l_X(i2)) + abs(s_l_X(i2) - s_l_X(i3)) + abs(s_l_X(i3) - s_l_X(i1)); 
    if sum_abs_x > 1
        index_row = [index_row; i]; 
    end
end
DT_s_l(index_row, :) = []; 

index_row = [];
for i = 1 : length(DT_h_s(:, 1))
    i1 = DT_h_s(i, 1); 
    i2 = DT_h_s(i, 2); 
    i3 = DT_h_s(i, 3); 
    C1 = [h_s_X(i1) h_s_Y(i1) h_s_Z(i1)/100];
    C2 = [h_s_X(i2) h_s_Y(i2) h_s_Z(i2)/100];   
    C3 = [h_s_X(i3) h_s_Y(i3) h_s_Z(i3)/100];
%     sum_abs_dist = abs(h_s_Y(i1) - h_s_Y(i2)) + abs(h_s_Y(i2) - h_s_Y(i3)) + abs(h_s_Y(i3) - h_s_Y(i1));
    sum_abs_dist = sqrt(sum((C1-C2).^2))+sqrt(sum((C2-C3).^2))+sqrt(sum((C3-C1).^2));
    if sum_abs_dist > 2.25
        index_row = [index_row; i]; 
    end
end
DT_h_s(index_row, :) = []; 
%     com_str = "trisurf(DT_s"+num2str(1)+",s_X(index_s"+num2str(1)+"),s_Y(index_s"+num2str(1)+"),sps(index_s"+num2str(1)+...
%         "),'FaceColor',colors(2,:),'FaceAlpha',0.2,'EdgeColor',colors(2,:));";

% eval(com_str);
%     triplot(DT_s1,s_X,s_Y,'Color',colors(2,:));
p4=trisurf(DT_ss, s_X, s_Y, sps, 'FaceColor', colors(2, :), 'FaceAlpha', 0.5, 'EdgeColor', 'none'); 
hold on
% trisurf(DT_s1,s_X,s_Y,spl_true,'FaceColor',colors(1,:),'FaceAlpha',0.2,'EdgeColor','none');
p3=trisurf(DT_sl, l_X, l_Y, spl, 'FaceColor', colors(1, :), 'FaceAlpha', 0.5, 'EdgeColor', 'none'); 
trisurf(DT_ss, h_X, h_Y, sph, 'FaceColor', colors(1, :), 'FaceAlpha', 0.5, 'EdgeColor', 'none'); 
trisurf(DT_l_h, l_h_X,l_h_Y, l_h_Z, 'FaceColor', colors(1, :), 'FaceAlpha', 0.5, 'EdgeColor', 'none'); 
trisurf(DT_s_l, s_l_X,s_l_Y, s_l_Z, 'FaceColor', colors(2, :), 'FaceAlpha', 0.5, 'EdgeColor', 'none'); 
trisurf(DT_h_s, h_s_X,h_s_Y, h_s_Z, 'FaceColor', colors(2, :), 'FaceAlpha', 0.5, 'EdgeColor', 'none'); 
grid on
hold on
axis([0 7 0 5.5 0 800]);
axis square

p1 = plot3(Y7(:, 4), ut_ZT7', Y7(:, 2), LineWidth = 1.5, Color = '#1F77B4'); 
p2 = plot3(Y19(:, 4), ut_ZT19', Y19(:, 2), LineWidth = 1.5, Color = '#D62728'); 
arrowPlot3D(Y7(:, 4), ut_ZT7', Y7(:, 2), 'numArrows', 2, 'color', '#1F77B4', 'headSize', 0.1, 'axis',[0 7 0 5.5 0 800],'LineVisible',0,...
    'SampleIndices', [500,1500]);
arrowPlot3D(Y19(:, 4), ut_ZT19', Y19(:, 2), 'numArrows', 2, 'color', '#D62728', 'headSize', 0.1, 'axis',[0 7 0 5.5 0 800],'LineVisible',0,...
    'SampleIndices', [600,1400]);
xlabel('DC_{PT} (10^3/\mul)', FontSize = 13); 
% ylabel('u(10^3\cdot\mul^{-1}\cdotd^{-1})',FontSize=13);
ylabel('u (10^3/(\mul·day))', FontSize = 13); 
zlabel('sCCL21 (pg/ml)', FontSize = 13); 
    set(gca, 'FontName', 'Arial'); 
    ax = gca; 
    ax.LineWidth = 1; % 加粗边框
    ax.TickLength = [0.008, 0.008]; % 前者设置y轴 tick的长度，后者设置x轴
    legend([p3, p4, p1, p2], {'Stable Steady States', 'Saddle Points', 'ZT7', 'ZT19'}, 'Location','best'); 



% DT_saddle = delaunay(Y(index2),sps);
% b_Other = boundary([Y(index1);Y(index3)],[spl;sph],0.1);
% YY_Other = [Y(index1);Y(index3)];
% ZZ_Other = [spl;sph];
% plot(YY_Other(b_Other),ZZ_Other(b_Other))
% b_Other1 = [b_Other(1:end);b_Other(1)];
% C_Other = [b_Other,b_Other1];
% DT_Other = delaunayTriangulation([Y(index1);Y(index3)],[spl;sph],C_Other);
% DT_Other = delaunay([Y(index1);Y(index3)],[spl;sph]);
% shp = alphaShape([Y(index1);Y(index3)],[spl;sph],10);
% plot(shp);
% tri = alphaTriangulation(shp);
% axis([0 7 0 750]);
% axis square
% nexttile(3)
% triplot(tri,[Y(index1);Y(index3)],[spl;sph]);
% triplot(DT,YY_s,ZZ_s);
% ConnectivityList
% trisurf(tri,[X(index1);X(index3)],[Y(index1);Y(index3)],[spl;sph],'FaceColor',colors(1,:),'FaceAlpha',0.3,'EdgeColor','none');
% p=trisurf(DT,XX_s,YY_s,ZZ_s,'FaceColor',colors(1,:),'FaceAlpha',0.2,'EdgeColor','none');
hold on
% trisurf(DT_saddle,X(index2),Y(index2),sps,'FaceColor',colors(2,:),'FaceAlpha',0.6,'EdgeColor','none');
hold off

nexttile(3)
H_state_index = find(Z(1, :, 3)); 
S_state_index = find(Z(1, :, 2)); 
L_state_index = find(~Z(1, :, 1)); 
H = Z(1, H_state_index, 3); 
S = Z(1, S_state_index, 2); 
L = Z(1, L_state_index, 1); 
Z_modify = (Z(1, min(S_state_index), 2) + Z(1, min(H_state_index), 3)) / 2; 
H = [Z_modify, Z(1, H_state_index, 3)]; 
S = [Z_modify, Z(1, S_state_index, 2)]; 
H_state_index = [min(H_state_index) - 1, H_state_index]; 
S_state_index = [min(S_state_index) - 1, S_state_index]; 

p1 = plot(DCPT(L_state_index), L, '-', 'Color', colors(1, :), 'LineWidth', 2); 
hold on
p2 = plot(DCPT(S_state_index), S, '--', 'Color', colors(2, :), 'LineWidth', 2); 
p3 = plot(DCPT(H_state_index), H, '-', 'Color', colors(1, :), 'LineWidth', 2); 
% p1=semilogy([6.89655;DC_PT_v;0.45], [67;fix_point_ensemble_DC(1:37,2);2.9]+0.1, '-', 'Color', colors(1,:), 'LineWidth', 1.5);
% hold on
% p2=semilogy([6.89655;DC_PT_v;0.45], [0.65;fix_point_ensemble_DC(1:37,1);2.9]+0.1, '--', 'Color', colors(2,:), 'LineWidth', 1.5);
% TEND = 50:50:4000;
% DC_PT_v = Y19(TEND,4);
% p3=semilogy([6.89655;DC_PT_v], [0;fix_point_ensemble_DC(:,3)]+0.1, '-', 'Color', colors(1,:), 'LineWidth', 1.5);
% DC_PT_v =Y19(:,4);
% p4=semilogy(DC_PT_v, Y19(:,3)+0.1, '-', 'Color', '#D62728', 'LineWidth', 1.5);
% DC_PT_v =Y7(:,4);
% p5=semilogy(DC_PT_v, Y7(:,3)+0.1, '-', 'Color', '#1F77B4', 'LineWidth', 1.5);
    set(gca, 'FontName', 'Arial'); 
    ax = gca; 
    ax.LineWidth = 1; % 加粗边框
    ax.TickLength = [0.008, 0.008]; % 前者设置y轴 tick的长度，后者设置x轴
    set(ax, 'TickDir', 'out')
    t = xlabel('DC_{PT} (10^3/\mul)', FontSize = 13); 
    t_position = t.Position; 
%     t.Position = [t_position(1),t_position(2)-0.5,t_position(3)];
    t = ylabel({'sCCL21 (pg/ml)'}, FontSize = 13); 
   
    title('u = 0', FontSize = 13); 
    legend([p1, p2], {'Stable Steady States', 'Saddle Points'},'Location','east'); 
    axis([0 7 -30 800]);
    axis square
    box off
    hold off



nexttile(4)
DCPT_i = 51; 
H_state_index = find(Z(:, DCPT_i, 3)); 
H_state_index1 = find(Z(:, DCPT_i, 1) > 200); 
S_state_index = find(Z(:, DCPT_i, 2)); 
L_state_index = find(Z(:, DCPT_i, 1) < 200); 

H = [Z(H_state_index, DCPT_i, 3); Z(H_state_index1, DCPT_i, 1)]; 
H_state_index = [H_state_index; H_state_index1]; 
S = Z(S_state_index, DCPT_i, 2); 
L = Z(L_state_index, DCPT_i, 1); 

p1 = plot(UT(L_state_index), L, '-', 'Color', colors(1, :), 'LineWidth', 2); 
hold on
p2 = plot(UT(S_state_index), S, '--', 'Color', colors(2, :), 'LineWidth', 2); 
p3 = plot(UT(H_state_index), H, '-', 'Color', colors(1, :), 'LineWidth', 2); 
    set(gca, 'FontName', 'Arial'); 
    ax = gca; 
    ax.LineWidth = 1; % 加粗边框
    ax.TickLength = [0.008, 0.008]; % 前者设置y轴 tick的长度，后者设置x轴
    set(ax, 'TickDir', 'out')
    t = xlabel('u (10^3/(\mul·day))', FontSize = 13); 
    t_position = t.Position; 
%     t.Position = [t_position(1),t_position(2)-0.5,t_position(3)];
    t = ylabel({'sCCL21(pg/ml)'}, FontSize = 13); 
   
    title('DC_{PT} = 0.5', FontSize = 13); 
    legend([p1, p2], {'Stable Steady States', 'Saddle Points'}, 'Location','east'); 
    axis([0 7 0 800]);
    axis square
    box off
    hold off
TT.TileSpacing = 'compact'; 
TT.Padding = 'compact'; 

end
function dydt = ode_DCCCL_new_ut(t, y, Par, Time_circadian_clock, ZT_initial)
circa_t_index = mod(floor((t * 24 + ZT_initial) / 24 * 1000), 1000) + 1; 
CLOCK_BMAL1_t = Time_circadian_clock(circa_t_index, 7); 
CCL21 = y(1); 
CCL19 = y(2); 
DC_LN = y(3); 
DC_PT = y(4); 
gamma_PTLN = Par(1); 
c_PTLN_cell = Par(2); 
K_CCL19_DC = Par(3); 
K_CCL21_DC = Par(4); 
alpha = Par(5); 
d_DC = Par(6); 
a_0_CCL = Par(7); 
a_DC_CCL = Par(8); 
K_DC_CCL = Par(9); 
d_CCL = Par(10); 

n_D = 2; 
n_C = 2; 

dydt = zeros(5, 1); 
% surface-immbolized CCL21
dydt(1) = a_0_CCL * CLOCK_BMAL1_t ^ 3 / (CLOCK_BMAL1_t ^ 3 + 0.1649 ^ 3) - d_CCL * CCL21; % x-CCL
% soluble CCL21 and CCL19
dydt(2) = a_DC_CCL * DC_LN ^ n_D / (DC_LN ^ n_D + K_DC_CCL ^ n_D) - d_CCL * CCL19; % x-CCL
% DC in LN
dydt(3) = gamma_PTLN * c_PTLN_cell * DC_PT * ( CCL21 ^ n_C / (K_CCL21_DC ^ n_C + CCL21 ^ n_C) + alpha * CCL19 ^ n_C / (K_CCL19_DC ^ n_C + CCL19 ^ n_C)) - d_DC * DC_LN; %y-DC
% DC in tissue
dydt(4) = - c_PTLN_cell * DC_PT * ( CCL21 ^ n_C / (K_CCL21_DC ^ n_C + CCL21 ^ n_C) + alpha * CCL19 ^ n_C / (K_CCL19_DC ^ n_C + CCL19 ^ n_C)) - d_DC * (DC_PT); 
dydt(5) = gamma_PTLN * c_PTLN_cell * DC_PT * CCL21 ^ n_C / (K_CCL21_DC ^ n_C + CCL21 ^ n_C); 
dydt(6) = gamma_PTLN * c_PTLN_cell * DC_PT * CCL21 ^ n_C / (K_CCL21_DC ^ n_C + CCL21 ^ n_C) - d_DC * y(6); 
dydt(7) = gamma_PTLN * c_PTLN_cell * DC_PT * alpha * CCL19 ^ n_C / (K_CCL19_DC ^ n_C + CCL19 ^ n_C) - d_DC * y(7); 
end
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



params = parseParams(varargin, defaults);

if size(X,1)>1  
    X = X';
    Y = Y';
    Z = Z';
end
V = [X;Y;Z];

if params.LineVisible==1

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




% total_length = sum(sqrt(diff(X).^2 + diff(Y).^2 + diff(Z).^2));
AMAX = [params.axis(2), params.axis(4), params.axis(6)];
AMIN = [params.axis(1), params.axis(3), params.axis(5)];
V_n = diag(1./(AMAX - AMIN)) * (V - diag(AMIN) * ones(3,length(X)));
X_n = V_n(1,:);
Y_n = V_n(2,:);
Z_n = V_n(3,:);




if params.SampleIndices == -1
    sample_indices = round(linspace(1, length(X), params.numArrows+2));
    sample_indices = sample_indices(2:end-1); 
else
    sample_indices = params.SampleIndices;
end


[arrowX, arrowY, arrowZ] = create3DArrow(params.headSize*params.scale);
DT = delaunay(arrowY, arrowZ);

for idx = sample_indices
    
    pos = [X_n(idx), Y_n(idx), Z_n(idx)];
    
    
    if idx == length(X)
        tangent = [X_n(idx)-X_n(idx-1), Y_n(idx)-Y_n(idx-1), Z_n(idx)-Z_n(idx-1)];
    else
        tangent = [X_n(idx+1)-X_n(idx), Y_n(idx+1)-Y_n(idx), Z_n(idx+1)-Z_n(idx)];
    end
    tangent = tangent / norm(tangent);
    
    
    R = vecRotation([1 0 0], tangent); 
    
   
    rotatedArrow = (R*[arrowX; arrowY; arrowZ]);
    pos = diag(pos)* ones(3, length(arrowX(1,:)));
    finalArrow = rotatedArrow + pos;
    finalArrow = diag(AMIN) * ones(3,length(arrowX(1,:))) + diag((AMAX - AMIN))*finalArrow;

    
   
%     fill3(finalArrow(1,:), finalArrow(2,:), finalArrow(3,:), ...
%           params.color, 'EdgeColor', 'none');

    trisurf(DT, finalArrow(1,:), finalArrow(2,:), finalArrow(3,:), 'FaceColor', params.color, 'FaceAlpha', 1, 'EdgeColor', 'none');

end

% axis vis3d;
view(3);
end


function params = parseParams(inputs, defaults)
params = defaults;
for i = 1:2:length(inputs)
    if isfield(defaults, inputs{i})
        params.(inputs{i}) = inputs{i+1};
    end
end
end


function [X, Y, Z] = create3DArrow(headSize)

theta = linspace(0, 2*pi, 64);
r = headSize/2;
z = r/2 * cos(theta);
y = r/2 * sin(theta);
x_base = zeros(size(z)) - r;


% X = [x_base, -headSize*1.5*ones(1,8), -headSize*1.5*ones(1,8)]';
% Y = [y, y*0.8, y*0.8]';
% Z = [z, z*0.8, z*0.8]';
X = x_base';
Y = y';
Z = z';

X = [X; 0]';
Y = [Y; 0]';
Z = [Z; 0]';

end


function R = vecRotation(vec1, vec2)

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


