function demo_fig4_fig5_fig6
Par_decay = xlsread('decay_list.xlsx');
Par_produce = xlsread('produce_list.xlsx');
Par_transition = xlsread('transition_list.xlsx');
Par_halfmax = xlsread('half_maximum_list.xlsx');
Par_Others = xlsread('Others_list.xlsx');
Par_promotion = xlsread('promotion_list.xlsx');
Par_clock = load('par_clock.csv');
Par_Clock_immune = xlsread('circa_immune_list.xlsx');
Time_circadian_clock = xlsread('Timing.xlsx');
load Par_ifam_ada_index.mat Par_iflam Par_ada
colors = lines(100);

dt = 0.001;
y0 = [0	0	0	0	0	0	0 ...
    0	0	0	0	0	0.0769532967279601	0	0	5.68586660722844e-06...
    0	0	0	0	0.995766071210120	0	0 ...
    41.2657516450204	0	0	400000.228733991	0	47.5877748302043	0 ...
    0	42.4288746428268	0	0	0	80.0544975496834...
    0	0	0	47.5877748302043	1.40067660106892	0 ...
    5.23905059949125e-06	0	0	0	0	0.128464137174231...
    0	0	0.391999999863231	0.958930572384777	0.983203150903828	1.88412743981574	0	0	0 ...
    0	0	0	0	3.83970865379364e-07	0	0 ...
    0	0	0.0133519230279894	0	0	0 ...
    0	0	1.54419834273011	1.33873453936487	1.62581918964826	2.36653249724990	2.30933994708977...
    10.2458891775095	0.518329265711375	0.00217250633691438 0.506885599518817	7.24910600997961	2.13186264381659	0.159473513507005];

% Endotoxin adjuvanted vaccine
% Par_Others(35)=1;

T7_0 = 0 : dt : 7 / 24+1;
T19_0 = 0 : dt : (19/ 24+1);
T7 = 7/24+1+dt : dt : 14;
T19 = (19/24+1+dt) : dt : 14;
% Tun_0 = 0:dt:14;

[~, Y7_0] = ode45(@(t, y)ODE_zero_circa_immune_response_fit_acute_version_T(t, y, Par_decay, Par_produce, ...
    Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,0), T7_0, y0);
[~, Y19_0] = ode45(@(t, y)ODE_zero_circa_immune_response_fit_acute_version_T(t, y, Par_decay, Par_produce, ...
    Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,0), T19_0, y0);
% % [~, Yun_0] = ode45(@(t, y)ODE_zero_circa_immune_response_fit_acute_version_T(t, y, Par_decay, Par_produce, ...
% %     Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,0), Tun_0, y0);
%
y7_0 = Y7_0(end, :);
y19_0 = Y19_0(end, :);
% LPS
p=4;
y7_0(1)=p;
y19_0(1)=p;
% CpG
% p=5;
% y7_0(1)=p;
% y19_0(1)=p;
% BMDC
% DC0 = 10^6/145000;
% y7_0(3)=DC0;
% y19_0(3)=DC0;
[~, Y7] = ode45(@(t, y)ODE_zero_circa_immune_response_fit_acute_version_T(t, y, Par_decay, Par_produce, ...
    Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,0), T7, y7_0);
[~, Y19] = ode45(@(t, y)ODE_zero_circa_immune_response_fit_acute_version_T(t, y, Par_decay, Par_produce, ...
    Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,0), T19, y19_0);
Y7=[Y7_0;Y7];
Y19=[Y19_0;Y19];
Tend=14;
T = 0:dt:Tend;
T = T(1:14000)*24;

inflam_ZT7 = inflammation(Y7(:,3),Y7(:,4),Y7(:,66),Y7(:,61),Y7(:,67),Y7(:,69),Par_iflam);
inflam_ZT19 = inflammation(Y19(:,3),Y19(:,4),Y19(:,66),Y19(:,61),Y19(:,67),Y19(:,69),Par_iflam);
ada_ZT7=Ada_response(Y7(:,34),Y7(:,35),Y7(:,68),Y7(:,39),Y7(:,64),Y7(:,65),Par_ada);
ada_ZT19=Ada_response(Y19(:,34),Y19(:,35),Y19(:,68),Y19(:,39),Y19(:,64),Y19(:,65),Par_ada);

i_index = [79,81,84,1,4,15,3,23,29,24,25,68];
y_vector = zeros(length(i_index),4);
axis_vector=zeros(length(i_index),4);
% subplot_index = [1,3,5,7,9,11,2,4,6,8,10,12];
subplot_index = [5,7,2,1,3,9,11,6,4,8,10,12];
index_turnover =1;
f=figure('Name','mechanism_LPS','NumberTitle','off');
f.Color = [1 1 1];
set(gcf, 'Position', [0 0 1300 800]); % 设置图形大小[宽度, 高度]
for i = 1:12
    Tend=7;
    if i ==3
        y_vector(i,:)=[0.1 0.1 0.225 0.225];
        axis_vector(i,:)=[0 Tend*24 0.1 0.225];
    elseif i==1
        y_vector(i,:)=[0.2 0.2 0.6 0.6];
        axis_vector(i,:)=[0 Tend*24 0.2 0.6];
    elseif i==2
        y_vector(i,:)=[0.3 0.3 0.6 0.6];
        axis_vector(i,:)=[0 Tend*24 0.3 0.6];
    elseif i==10
        maximum=1.2*log10(max(Y7(:,26))+1);
        y_vector(i,:)=[0 0 maximum+1 maximum+1];
        axis_vector(i,:)=[0 Tend*24 0 maximum+1];
    elseif i==11
        maximum=1.2*max(Y7(:,37));
        y_vector(i,:)=[0 0 maximum+1 maximum+1];
        axis_vector(i,:)=[0 Tend*24 0 maximum+1];
    elseif i==7
%         maximum=1.2*max(Y7(:,3));
        if index_turnover == 1
            maximum=1.2*max(Y7(:,3));
        else maximum=1.2*max(inflam_ZT7);
        end
        y_vector(i,:)=[0 0 maximum maximum];
        axis_vector(i,:)=[0 Tend*24 0 maximum];
    elseif i==12
        Tend=14;
        if index_turnover == 1
        maximum=1.2*max(Y7(:,68));
        else
        maximum=40;
        end
        y_vector(i,:)=[0 0 maximum maximum];
        axis_vector(i,:)=[0 Tend*24 0 maximum];

    elseif i==4
        if Par_Others(42)==0
            maximum=1.2*max(Y7(:,i_index(i)));
            y_vector(i,:)=[0 0 maximum+1 maximum+1];
            axis_vector(i,:)=[0 Tend*24 0 maximum+1];
        else
            maximum=1.2*max(Y7(:,i_index(i)));
            y_vector(i,:)=[0 0 maximum+1 maximum+1];
            axis_vector(i,:)=[0 Tend*24 0 maximum+1];
        end
   elseif i==8
        if Par_Others(42)==0
            maximum=1.2*max(Y7(:,i_index(i)));
            y_vector(i,:)=[0 0 maximum+1 maximum+1];
            axis_vector(i,:)=[0 Tend*24 0 maximum+1];
        else
            maximum=1.2*max(Y7(:,i_index(i)));
            y_vector(i,:)=[0 0 maximum+1 maximum+1];
            axis_vector(i,:)=[0 Tend*24 0 maximum+1];
        end
    else
        maximum=1.2*max(Y7(:,i_index(i)));
        y_vector(i,:)=[0 0 maximum+1 maximum+1];
        axis_vector(i,:)=[0 Tend*24 0 maximum+1];
    end
end
% if index_turnover ==1
title_cell = {'CRYs','RORs','CLOCK:BMAL1','Endotoxin-adjuvanted Vaccine (Ag)','Macrophages in Tissues (M\Phi_{PT})','TNF-\alpha in Tissues (TNF-\alpha_{PT})',...
    'DCs in Tissues (DC_{PT})','DCs in dLNs (DC_{LN})','CCL21','Helper T Cells (Th)','Germinal Center B Cells (B_{GC})','Antibodies in Blood (Ab)'};
Ylabel_cell = {'(nM)','(nM)','(nM)','(mg/kg)','(10^3/\mul)','(pg/ml)',...
    '(10^3/\mul)','(10^3/\mul)','(pg/ml)','(10^3/\mul)','(10^3/\mul)','(\mug/ml)'};
eps = 0.0001;
ylim_cell = {[0.2 0.6+eps],[0.3 0.6+eps],[0.1 0.25+eps], [0 6+eps],[0 400+eps],[0 60000+eps],[0 8+eps],[0 6+eps],...
    [0 600+eps],[0 1000+eps],[0 100],[0 30+eps]};
y_vector_highest = [0.6, 0.6, 0.25, 6, 400, 60000, 8, 6, 600, 1000, 100, 40];
y_vector(:,3) = y_vector_highest;
y_vector(:,4) = y_vector_highest;
yticks_cell = {[0.2 0.4 0.6],[0.3, 0.4, 0.5, 0.6],[0.1,0.15,0.2, 0.25],[0 2 4 6],[0,200,400],[0 20000 40000 60000],[0 4 8],...
    [0,2,4,6],[0,200,400,600],[0 500 1000], [0 50 100], [0 10 20 30]};
% end

for i = 1:12
    x_i = mod(subplot_index(i)-1,2)+1;
    y_i = ceil(subplot_index(i)/2);
    al = 0.42;
    bl = 0.09;
    ax=subplot(6,2,subplot_index(i));
    ax.Position=[0.25+(x_i-1)*0.5, 1-(0.12+0.155*(y_i-1)),0,0]+[-al/2,-bl/2,0,0]+[0.02,0.03,0,0]+[0,0,al,bl];

    if i == 10
    plot(T(1:7000), Y7(1:7000, 26), '-', 'Color', colors(1, :), 'LineWidth', 1.5);
    hold on
    plot(T(1:7000), Y19(1:7000, 26), '-', 'Color', colors(2, :), 'LineWidth', 1.5);
    symlog_2('y',0,'Ylim',[0 1000]);


    elseif i ==11
    plot(T(1:7000), Y7(1:7000, 37), '-', 'Color', colors(1, :), 'LineWidth', 1.5);
    hold on
    plot(T(1:7000), Y19(1:7000, 37), '-', 'Color', colors(2, :), 'LineWidth', 1.5);
    elseif i==9
    plot(T(1:7000), Y7(1:7000,i_index(i)), '-', 'Color', colors(1, :), 'LineWidth', 1.5);
    hold on
    plot(T(1:7000), Y19(1:7000,i_index(i)), '-', 'Color', colors(2, :), 'LineWidth', 1.5);
%     plot(Tun_0(1:7000)*24, Yun_0(1:7000,i_index(i)), '--', 'Color', colors(2, :), 'LineWidth', 1.5);
    elseif i==7
        if index_turnover==1
    plot(T(1:7000), Y7(1:7000,i_index(i)), '-', 'Color', colors(1, :), 'LineWidth', 1.5);
    hold on
    plot(T(1:7000), Y19(1:7000,i_index(i)), '-', 'Color', colors(2, :), 'LineWidth', 1.5);
        else
    plot(T(1:7000), inflam_ZT7(1:7000), '-', 'Color', colors(1, :), 'LineWidth', 1.5);
    hold on
    plot(T(1:7000), inflam_ZT19(1:7000), '-', 'Color', colors(2, :), 'LineWidth', 1.5);
        end
    elseif i==4
        if Par_Others(42)==0
           plot(T(1:7000), Y7(1:7000,i_index(i)), '-', 'Color', colors(1, :), 'LineWidth', 1.5);
           hold on
           plot(T(1:7000), Y19(1:7000,i_index(i)), '-', 'Color', colors(2, :), 'LineWidth', 1.5);
        else
           plot(T(1:7000), Y7(1:7000,i_index(i)), '-', 'Color', colors(1, :), 'LineWidth', 1.5);
           hold on
           plot(T(1:7000), Y19(1:7000,i_index(i)), '-', 'Color', colors(2, :), 'LineWidth', 1.5);
        end
    elseif i==12
        if index_turnover==1
    plot(T, Y7(:,i_index(i)), '-', 'Color', colors(1, :), 'LineWidth', 1.5);
    hold on
    plot(T, Y19(:,i_index(i)), '-', 'Color', colors(2, :), 'LineWidth', 1.5);
        else
    plot(T, ada_ZT7, '-', 'Color', colors(1, :), 'LineWidth', 1.5);
    hold on
    plot(T, ada_ZT19, '-', 'Color', colors(2, :), 'LineWidth', 1.5);
        end
    else
    p114=plot(T(1:7000), Y7(1:7000,i_index(i)), '-', 'Color', colors(1, :), 'LineWidth', 1.5);
    hold on
    p514=plot(T(1:7000), Y19(1:7000,i_index(i)), '-', 'Color', colors(2, :), 'LineWidth', 1.5);
    end
    if i==12
        for j=1:14
        p_command_string_1 = "p"+num2str(i)+"_"+num2str(j)+"=patch([24*"+num2str(j)+"-12,24*"+num2str(j)+",24*"+num2str(j)+...
            ", 24*"+num2str(j)+"-12],y_vector("+num2str(i)+",:),'k');";
        p_command_string_2 = "p"+num2str(i)+"_"+num2str(j)+".FaceColor =[0.5 0.5 0.5];";
        p_command_string_3 = "p"+num2str(i)+"_"+num2str(j)+".FaceVertexAlphaData =0.5;";
        p_command_string_4 = "p"+num2str(i)+"_"+num2str(j)+".FaceAlpha = 'flat' ;";
        p_command_string_5 = "p"+num2str(i)+"_"+num2str(j)+".EdgeColor = 'none';";
        eval(p_command_string_1);
        eval(p_command_string_2);
        eval(p_command_string_3);
        eval(p_command_string_4);
        eval(p_command_string_5);
        end
    else
        for j = 1:7
           
        p_command_string_1 = "p"+num2str(i)+"_"+num2str(j)+"=patch([24*"+num2str(j)+"-12,24*"+num2str(j)+",24*"+num2str(j)+...
            ", 24*"+num2str(j)+"-12],y_vector("+num2str(i)+",:),'k');";
        p_command_string_2 = "p"+num2str(i)+"_"+num2str(j)+".FaceColor =[0.5 0.5 0.5];";
        p_command_string_3 = "p"+num2str(i)+"_"+num2str(j)+".FaceVertexAlphaData =0.5;";
        p_command_string_4 = "p"+num2str(i)+"_"+num2str(j)+".FaceAlpha = 'flat' ;";
        p_command_string_5 = "p"+num2str(i)+"_"+num2str(j)+".EdgeColor = 'none';";
        eval(p_command_string_1);
        eval(p_command_string_2);
        eval(p_command_string_3);
        eval(p_command_string_4);
        eval(p_command_string_5);
        end
    end
    if index_turnover==1
    if i~=10
    axis(axis_vector(i,:));
    ylim(ylim_cell{i});
    yticks(yticks_cell{i});
    else
        xlim([0 24*7])
    end
    end
    if index_turnover==0
        if i==7
            lim = [0 4+eps];
            ytic = [0 2 4];
            axis(axis_vector(i,:));
            ylim(lim);
            yticks(ytic);
        end
        if i==12
            lim = [0 40+eps];
            ytic = [0 20 40];
            axis(axis_vector(i,:));
            ylim(lim);
            yticks(ytic);
        end
    end
    set(gca, 'FontName', 'Arial');
    ax = gca;
    ax.LineWidth = 1.5; 
    ax.TickLength = [0.004, 0.004]; 
    ax.FontSize = 14;
    set(ax, 'TickDir', 'out');
    box off;
    hold off;
    if i ==1
        hl=legend([p114,p514],{'ZT7    ','ZT19'},'Orientation','horizontal',FontSize=14);
        hl.Box = 'off';
    end
    tit=title(title_cell{i});
    tit.FontSize = 16;
    t=ylabel(Ylabel_cell{i});
    t_position=t.Position;
    t.Position = [-11,t_position(2),t_position(3)];
    if i ==12
       t.Position = [-18,t_position(2),t_position(3)];
    else
       t.Position = [-9,t_position(2),t_position(3)];
    end
    t.FontSize = 16;
    if ismember(i,[7,12])
        xl=xlabel('Time (hour)');
        xl.FontSize = 16;
        xl.Position=[xl.Position(1) xl.Position(2)*0.6 xl.Position(3)];
    end

    if i<12
    xticks(0:12:24*Tend);
    xtickangle(0);
%     zero_str = '0';
%     tw_str = '12';
%     ZT_ticklabel={zero_str,tw_str};
%     xticklabel_str = {'ZT0','12'};
    else
    xticks(0:24:24*Tend);
    xtickangle(0);
%     zero_str = '0';

%     ZT_ticklabel={zero_str};
%     xticklabel_str = {'ZT0'};
    end

    if ismember(i,[1,2,3,4,5,6,8,9,10])
        xticklabels({});
    end
%     for k = 1:Tend-1
%         if k==Tend-1
% %         xticklabel_str=[xticklabel_str,ZT_ticklabel,{'0'}];
%         else
% %         xticklabel_str=[xticklabel_str,ZT_ticklabel];
%         end
%     end
%     xticklabels(xticklabel_str);
end

% Nucleic-acid adjuvanted vaccine
% endotoxin-Nucleic-acid controller
Par_Others(35)=1;

T7_0 = 0 : dt : 7 / 24+1;
T19_0 = 0 : dt : (19/ 24+1);
T7 = 7/24+1+dt : dt : 14;
T19 = (19/24+1+dt) : dt : 14;
% Tun_0 = 0:dt:14;

[~, Y7_0] = ode45(@(t, y)ODE_zero_circa_immune_response_fit_acute_version_T(t, y, Par_decay, Par_produce, ...
    Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,0), T7_0, y0);
[~, Y19_0] = ode45(@(t, y)ODE_zero_circa_immune_response_fit_acute_version_T(t, y, Par_decay, Par_produce, ...
    Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,0), T19_0, y0);
% % [~, Yun_0] = ode45(@(t, y)ODE_zero_circa_immune_response_fit_acute_version_T(t, y, Par_decay, Par_produce, ...
% %     Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,0), Tun_0, y0);
%
y7_0 = Y7_0(end, :);
y19_0 = Y19_0(end, :);
% LPS
% p=4;
% y7_0(1)=p;
% y19_0(1)=p;
% CpG
p=5;
y7_0(1)=p;
y19_0(1)=p;
% BMDC
% DC0 = 10^6/145000;
% y7_0(3)=DC0;
% y19_0(3)=DC0;
[~, Y7] = ode45(@(t, y)ODE_zero_circa_immune_response_fit_acute_version_T(t, y, Par_decay, Par_produce, ...
    Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,0), T7, y7_0);
[~, Y19] = ode45(@(t, y)ODE_zero_circa_immune_response_fit_acute_version_T(t, y, Par_decay, Par_produce, ...
    Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,0), T19, y19_0);
Y7=[Y7_0;Y7];
Y19=[Y19_0;Y19];
Tend=14;
T = 0:dt:Tend;
T = T(1:14000)*24;

inflam_ZT7 = inflammation(Y7(:,3),Y7(:,4),Y7(:,66),Y7(:,61),Y7(:,67),Y7(:,69),Par_iflam);
inflam_ZT19 = inflammation(Y19(:,3),Y19(:,4),Y19(:,66),Y19(:,61),Y19(:,67),Y19(:,69),Par_iflam);
ada_ZT7=Ada_response(Y7(:,34),Y7(:,35),Y7(:,68),Y7(:,39),Y7(:,64),Y7(:,65),Par_ada);
ada_ZT19=Ada_response(Y19(:,34),Y19(:,35),Y19(:,68),Y19(:,39),Y19(:,64),Y19(:,65),Par_ada);

i_index = [1,84,3,29,23,26,37,68,99,100];
y_vector = zeros(length(i_index),4);
axis_vector=zeros(length(i_index),4);
subplot_index = [1,3,5,7,2,4,6,8,9,10];
f=figure('Name','mechanism_CpG','NumberTitle','off');
f.Color = [1 1 1];
set(gcf, 'Position', [0 0 1300 800]); % 设置图形大小[宽度, 高度]
% TTT = tiledlayout(5, 2);
% set(gcf,'unit','centimeters','position',[1,1,20,25]);
for i = 1:length(i_index)
    Tend=7;
    if i ==2
        y_vector(i,:)=[0.1 0.1 0.225 0.225];
        axis_vector(i,:)=[0 Tend*24 0.1 0.225];
    elseif i==8
        Tend=14;
        maximum=1.2*max(Y19(:,i_index(i)));
        y_vector(i,:)=[0 0 maximum+1 maximum+1];
        axis_vector(i,:)=[0 Tend*24 0 maximum+1];
    elseif i==6
        maximum=1.2*max(Y19(:,i_index(i)));
        y_vector(i,:)=[0 0 maximum maximum];
        axis_vector(i,:)=[0 Tend*24 0 maximum];
    elseif i==9
        maximum=1.2*max(inflam_ZT19);
        y_vector(i,:)=[0 0 maximum maximum];
        axis_vector(i,:)=[0 Tend*24 0 maximum];
    elseif i==10
        Tend=14;
        maximum=1.2*max(ada_ZT19);
        y_vector(i,:)=[0 0 maximum maximum];
        axis_vector(i,:)=[0 Tend*24 0 maximum];
    else
        maximum=1.2*max(Y19(:,i_index(i)));
        y_vector(i,:)=[0 0 maximum+1 maximum+1];
        axis_vector(i,:)=[0 Tend*24 0 maximum+1];
    end
end
title_cell = {'Nucleic-acid-adjuvanted Vaccine (Ag)','CLOCK:BMAL1','DCs in Tissues (DC_{PT})','CCL21',...
    'DCs in dLNs (DC_{LN})','Helper T Cells (Th)','Germinal Center B Cells (B_{GC})','Antibodies in Blood (Ab)','Inflammation','Adaptive Response'};
Ylabel_cell = {'(mg/kg)','(nM)','(10^3/\mul)','(pg/ml)',...
    '(10^3/\mul)','(10^3/\mul)','(10^3/\mul)','(\mug/ml)','',''};
eps = 0.0001;
ylim_cell = {[0 6+eps],[0.1 0.25+eps],[0 12+eps], [0 800+eps],[0 15+eps],[0 600+eps],[0 200+eps],[0 40+eps],...
    [0 3+eps],[0 60+eps]};
y_vector_highest = zeros(1,10);
for i = 1:10
y_vector_highest(i) = ylim_cell{i}(2);
end
y_vector(:,3) = y_vector_highest;
y_vector(:,4) = y_vector_highest;
% max_matrix = [4; 0.25; 10; 300; 2;  75; 75; 30;3;4];
%
yticks_cell = {[0 3 6],[0.1 0.15,0.2,0.25],[0 6 12],...
    0:400:800,[0 5 10 15],0:200:600,0:100:200,0:20:40,0:3,[0 20 60]};

for i = 1:length(i_index)
    x_i = mod(subplot_index(i)-1,2)+1;
    y_i = ceil(subplot_index(i)/2);
    al = 0.42;
    bl = 0.09;
    ax=subplot(6,2,subplot_index(i));
    ax.Position=[0.25+(x_i-1)*0.5, 1-(0.12+0.155*(y_i-1)),0,0]+[-al/2,-bl/2,0,0]+[0.02,0.03,0,0]+[0,0,al,bl];
%     nexttile(subplot_index(i));
    if i == 6
    plot(T(1:7000), Y7(1:7000, 26), '-', 'Color', colors(1, :), 'LineWidth', 1.5);
    hold on
    plot(T(1:7000), Y19(1:7000, 26), '-', 'Color', colors(2, :), 'LineWidth', 1.5);
%     symlog('y',-0.2);
    elseif i==8
    plot(T, Y7(:,68), '-', 'Color', colors(1, :), 'LineWidth', 1.5);
    hold on
    plot(T, Y19(:,68), '-', 'Color', colors(2, :), 'LineWidth', 1.5);
    elseif i==9
    plot(T(1:7000), inflam_ZT7(1:7000), '-', 'Color', colors(1, :), 'LineWidth', 1.5);
    hold on
    plot(T(1:7000), inflam_ZT19(1:7000), '-', 'Color', colors(2, :), 'LineWidth', 1.5);
    elseif i==10
    plot(T, ada_ZT7, '-', 'Color', colors(1, :), 'LineWidth', 1.5);
    hold on
    plot(T, ada_ZT19, '-', 'Color', colors(2, :), 'LineWidth', 1.5);

    else
    p114=plot(T(1:7000), Y7(1:7000,i_index(i)), '-', 'Color', colors(1, :), 'LineWidth', 1.5);
    hold on
    p514=plot(T(1:7000), Y19(1:7000,i_index(i)), '-', 'Color', colors(2, :), 'LineWidth', 1.5);
    end
    if ismember(i,[8,10])
        for j=1:14
        p_command_string_1 = "p"+num2str(i)+"_"+num2str(j)+"=patch([24*"+num2str(j)+"-12,24*"+num2str(j)+",24*"+num2str(j)+...
            ", 24*"+num2str(j)+"-12],y_vector("+num2str(i)+",:),'k');";
        p_command_string_2 = "p"+num2str(i)+"_"+num2str(j)+".FaceColor =[0.5 0.5 0.5];";
        p_command_string_3 = "p"+num2str(i)+"_"+num2str(j)+".FaceVertexAlphaData =0.5;";
        p_command_string_4 = "p"+num2str(i)+"_"+num2str(j)+".FaceAlpha = 'flat' ;";
        p_command_string_5 = "p"+num2str(i)+"_"+num2str(j)+".EdgeColor = 'none';";
        eval(p_command_string_1);
        eval(p_command_string_2);
        eval(p_command_string_3);
        eval(p_command_string_4);
        eval(p_command_string_5);
        end
    else
        for j = 1:7
        p_command_string_1 = "p"+num2str(i)+"_"+num2str(j)+"=patch([24*"+num2str(j)+"-12,24*"+num2str(j)+",24*"+num2str(j)+...
            ", 24*"+num2str(j)+"-12],y_vector("+num2str(i)+",:),'k');";
        p_command_string_2 = "p"+num2str(i)+"_"+num2str(j)+".FaceColor =[0.5 0.5 0.5];";
        p_command_string_3 = "p"+num2str(i)+"_"+num2str(j)+".FaceVertexAlphaData =0.5;";
        p_command_string_4 = "p"+num2str(i)+"_"+num2str(j)+".FaceAlpha = 'flat' ;";
        p_command_string_5 = "p"+num2str(i)+"_"+num2str(j)+".EdgeColor = 'none';";
        eval(p_command_string_1);
        eval(p_command_string_2);
        eval(p_command_string_3);
        eval(p_command_string_4);
        eval(p_command_string_5);
        end
    end
    % %   axis_vector(i,:)=[0 Tend*24 0 max_matrix(i)];
        axis(axis_vector(i,:));
        ylim(ylim_cell{i});
        yticks(yticks_cell{i});
    set(gca, 'FontName', 'Arial');
    ax = gca;
    ax.LineWidth = 1.5; 
    ax.TickLength = [0.004, 0.004];
    ax.FontSize = 14;
    set(ax, 'TickDir', 'out');
    box off;
    hold off;
    if i ==1
        hl=legend([p114,p514],{'ZT7','ZT19'},'Orientation','horizontal',FontSize=14);
        hl.Box = 'off';
    end
    tit=title(title_cell{i});
    tit.FontSize = 16;
    t=ylabel(Ylabel_cell{i});
    t_position=t.Position;
%     t.Position = [-6,t_position(2),t_position(3)];
    if i ==8
       t.Position = [-22,t_position(2),t_position(3)];
    else
       t.Position = [-11,t_position(2),t_position(3)];
    end
    t.FontSize = 16;
    if ismember(i,[9,10])
        xl=xlabel('Time (hour)');
        xl.FontSize = 16;
        xl.Position=[xl.Position(1) xl.Position(2)*0.3 xl.Position(3)];
    end
    if ismember(i,[8,10])
    xticks(0:24:24*Tend);
    xtickangle(0);
%     zero_str = '0';

%     ZT_ticklabel={zero_str};
%     xticklabel_str = {'ZT0'};
    else
    xticks(0:12:24*Tend);
    xtickangle(0);
%     zero_str = '0';
%     tw_str = '12';
%     ZT_ticklabel={zero_str,tw_str};
%     xticklabel_str = {'ZT0','12'};
    end
%     for k = 1:Tend-1
%         if k==Tend-1
% %         xticklabel_str=[xticklabel_str,ZT_ticklabel,{'0'}];
%         else
% %         xticklabel_str=[xticklabel_str,ZT_ticklabel];
%         end
%     end
%     xticklabels(xticklabel_str);
    if ismember(i,[1,2,3,4,5,6,8])
        xticklabels({});
    end
end

% DC vaccine

% Par_Others(35)=1;

T7_0 = 0 : dt : 7 / 24+1;
T19_0 = 0 : dt : (19/ 24+1);
T7 = 7/24+1+dt : dt : 14;
T19 = (19/24+1+dt) : dt : 14;
% Tun_0 = 0:dt:14;

[~, Y7_0] = ode45(@(t, y)ODE_zero_circa_immune_response_fit_acute_version_T(t, y, Par_decay, Par_produce, ...
    Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,0), T7_0, y0);
[~, Y19_0] = ode45(@(t, y)ODE_zero_circa_immune_response_fit_acute_version_T(t, y, Par_decay, Par_produce, ...
    Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,0), T19_0, y0);
% % [~, Yun_0] = ode45(@(t, y)ODE_zero_circa_immune_response_fit_acute_version_T(t, y, Par_decay, Par_produce, ...
% %     Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,0), Tun_0, y0);
%
y7_0 = Y7_0(end, :);
y19_0 = Y19_0(end, :);
% LPS
% p=4;
% y7_0(1)=p;
% y19_0(1)=p;
% CpG
% p=5;
% y7_0(1)=p;
% y19_0(1)=p;
% BMDC
DC0 = 10^6/145000;
y7_0(3)=DC0;
y19_0(3)=DC0;
[~, Y7] = ode45(@(t, y)ODE_zero_circa_immune_response_fit_acute_version_T(t, y, Par_decay, Par_produce, ...
    Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,0), T7, y7_0);
[~, Y19] = ode45(@(t, y)ODE_zero_circa_immune_response_fit_acute_version_T(t, y, Par_decay, Par_produce, ...
    Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,0), T19, y19_0);
Y7=[Y7_0;Y7];
Y19=[Y19_0;Y19];
Tend=14;
T = 0:dt:Tend;
T = T(1:14000)*24;

inflam_ZT7 = inflammation(Y7(:,3),Y7(:,4),Y7(:,66),Y7(:,61),Y7(:,67),Y7(:,69),Par_iflam);
inflam_ZT19 = inflammation(Y19(:,3),Y19(:,4),Y19(:,66),Y19(:,61),Y19(:,67),Y19(:,69),Par_iflam);
ada_ZT7=Ada_response(Y7(:,34),Y7(:,35),Y7(:,68),Y7(:,39),Y7(:,64),Y7(:,65),Par_ada);
ada_ZT19=Ada_response(Y19(:,34),Y19(:,35),Y19(:,68),Y19(:,39),Y19(:,64),Y19(:,65),Par_ada);

i_index = [79,81,84,1,4,15,3,23,29,24,25,68];
i_index = [3,84,29,23,26,37,68,100];
y_vector = zeros(length(i_index),4);
axis_vector=zeros(length(i_index),4);

subplot_index = [1,3,5,7,2,4,6,8];
f=figure('Name','mechanism_BMDC','NumberTitle','off');
f.Color = [1 1 1];
set(gcf, 'Position', [0 0 1300 800]); 
for i = 1:length(i_index)
    Tend=7;
    if i ==2
        y_vector(i,:)=[0.1 0.1 0.225 0.225];
        axis_vector(i,:)=[0 Tend*24 0.1 0.225];
    elseif i==7
        Tend=14;
        maximum=1.2*max(Y7(:,i_index(i)));
        y_vector(i,:)=[0 0 maximum+1 maximum+1];
        axis_vector(i,:)=[0 Tend*24 0 maximum+1];
    elseif i==5
        maximum=1.2*max(Y7(:,i_index(i)));
        y_vector(i,:)=[0 0 maximum maximum];
        axis_vector(i,:)=[0 Tend*24 0 maximum];
    elseif i==8
        Tend=14;
        maximum=1.2*max(ada_ZT7);
        y_vector(i,:)=[0 0 maximum maximum];
        axis_vector(i,:)=[0 Tend*24 0 maximum];
    else
        maximum=1.2*max(Y7(:,i_index(i)));
        y_vector(i,:)=[0 0 maximum+1 maximum+1];
        axis_vector(i,:)=[0 Tend*24 0 maximum+1];
    end
end

title_cell = {'DC Vaccine (BMDCs)','CLOCK:BMAL1','CCL21',...
    'DCs in dLNs (DC_{LN})','Helper T Cells (Th)','Germinal Center B cells (B_{GC})','Antibodies in Blood (Ab)','Adaptive Response'};
Ylabel_cell = {'(10^3/\mul)','(nM)','(pg/ml)',...
    '(10^3/\mul)','(10^3/\mul)','(10^3/\mul)','(\mug/ml)',''};
eps = 0.0001;
ylim_cell = {[0 10+eps],[0.1 0.25+eps],[0 750+eps], [0 10+eps],[0 1200+eps],[0 150+eps],[0 40+eps],[0 60+eps]};
y_vector_highest = [10; 0.25; 750; 10; 1200; 150; 50; 60]+eps;
y_vector(:,3) = y_vector_highest;
y_vector(:,4) = y_vector_highest;
max_matrix = [10,0.25,750,10,600,150,50,20];

yticks_cell = {[0 5 10],[0.1 0.15,0.2,0.25],0:250:750,...
    [0 5 10],0:400:1200,0:50:150,[0 20 40],[0 20 40 60]};
for i = 1:length(i_index)
    x_i = mod(subplot_index(i)-1,2)+1;
    y_i = ceil(subplot_index(i)/2);
    al = 0.42;
    bl = 0.09;
    ax=subplot(6,2,subplot_index(i));
    ax.Position=[0.25+(x_i-1)*0.5, 1-(0.12+0.155*(y_i-1)),0,0]+[-al/2,-bl/2,0,0]+[0.02,0.03,0,0]+[0,0,al,bl];

    if i == 5
    plot(T(1:7000), Y7(1:7000, 26), '-', 'Color', colors(1, :), 'LineWidth', 1.5);
    hold on
    plot(T(1:7000), Y19(1:7000, 26), '-', 'Color', colors(2, :), 'LineWidth', 1.5);
    ylim(ylim_cell{i});
    elseif i==7
    plot(T, Y7(:,68), '-', 'Color', colors(1, :), 'LineWidth', 1.5);
    hold on
    plot(T, Y19(:,68), '-', 'Color', colors(2, :), 'LineWidth', 1.5);
    ylim(ylim_cell{i});
    elseif i==8
    plot(T, ada_ZT7, '-', 'Color', colors(1, :), 'LineWidth', 1.5);
    hold on
    plot(T, ada_ZT19, '-', 'Color', colors(2, :), 'LineWidth', 1.5);
    ylim(ylim_cell{i});
    else
    p114=plot(T(1:7000), Y7(1:7000,i_index(i)), '-', 'Color', colors(1, :), 'LineWidth', 1.5);
    hold on
    p514=plot(T(1:7000), Y19(1:7000,i_index(i)), '-', 'Color', colors(2, :), 'LineWidth', 1.5);
    ylim(ylim_cell{i});
    end
    if ismember(i,[7,8])
        for j=1:14
        p_command_string_1 = "p"+num2str(i)+"_"+num2str(j)+"=patch([24*"+num2str(j)+"-12,24*"+num2str(j)+",24*"+num2str(j)+...
            ", 24*"+num2str(j)+"-12],y_vector("+num2str(i)+",:),'k');";
        p_command_string_2 = "p"+num2str(i)+"_"+num2str(j)+".FaceColor =[0.5 0.5 0.5];";
        p_command_string_3 = "p"+num2str(i)+"_"+num2str(j)+".FaceVertexAlphaData =0.5;";
        p_command_string_4 = "p"+num2str(i)+"_"+num2str(j)+".FaceAlpha = 'flat' ;";
        p_command_string_5 = "p"+num2str(i)+"_"+num2str(j)+".EdgeColor = 'none';";
        eval(p_command_string_1);
        eval(p_command_string_2);
        eval(p_command_string_3);
        eval(p_command_string_4);
        eval(p_command_string_5);
        end
    else
        for j = 1:7
        p_command_string_1 = "p"+num2str(i)+"_"+num2str(j)+"=patch([24*"+num2str(j)+"-12,24*"+num2str(j)+",24*"+num2str(j)+...
            ", 24*"+num2str(j)+"-12],y_vector("+num2str(i)+",:),'k');";
        p_command_string_2 = "p"+num2str(i)+"_"+num2str(j)+".FaceColor =[0.5 0.5 0.5];";
        p_command_string_3 = "p"+num2str(i)+"_"+num2str(j)+".FaceVertexAlphaData =0.5;";
        p_command_string_4 = "p"+num2str(i)+"_"+num2str(j)+".FaceAlpha = 'flat' ;";
        p_command_string_5 = "p"+num2str(i)+"_"+num2str(j)+".EdgeColor = 'none';";
        eval(p_command_string_1);
        eval(p_command_string_2);
        eval(p_command_string_3);
        eval(p_command_string_4);
        eval(p_command_string_5);
        end
    end

    axis(axis_vector(i,:));
    ylim(ylim_cell{i});
    yticks(yticks_cell{i});
    set(gca, 'FontName', 'Arial');
    ax = gca;
    ax.LineWidth = 1.5; % 
    ax.TickLength = [0.004, 0.004]; % 
    ax.FontSize = 14;
    set(ax, 'TickDir', 'out');
    box off;
    hold off;
    if i ==1
        hl=legend([p114,p514],{'ZT7','ZT19'},'Orientation','horizontal',FontSize=14);
        hl.Box = 'off';
    end
    tit=title(title_cell{i});
    tit.FontSize = 16;
    t=ylabel(Ylabel_cell{i});
    t_position=t.Position;
    t.Position = [-11,t_position(2),t_position(3)];
    if i ==7
       t.Position = [-18,t_position(2),t_position(3)];
    else
       t.Position = [-9,t_position(2),t_position(3)];
    end
    t.FontSize = 16;
    if ismember(i,[4,8])
        xl=xlabel('Time (hour)');
        xl.FontSize = 16;
        xl.Position=[xl.Position(1) xl.Position(2)*0.6 xl.Position(3)];
    end
    if ismember(i,[7,8])
    xticks(0:24:24*Tend);
    xtickangle(0);
%     zero_str = '0';
%
%     ZT_ticklabel={zero_str};
%     xticklabel_str = {'0'};
    else
    xticks(0:12:24*Tend);
    xtickangle(0);
%     zero_str = '0';
%     tw_str = '12';
%     ZT_ticklabel={zero_str,tw_str};
%     xticklabel_str = {'0','12'};
%     end
%     for k = 1:Tend-1
%         if k==Tend-1
%         xticklabel_str=[xticklabel_str,ZT_ticklabel,{'0'}];
%         else
%         xticklabel_str=[xticklabel_str,ZT_ticklabel];
%         end
%     end
%     xticklabels(xticklabel_str);
    end
    if ismember(i,[1,2,3,5,7])
        xticklabels({});
    end
end

end

function dydt = ODE_zero_circa_immune_response_fit_acute_version_T(t, y, Par_decay, Par_produce, ...
    Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune, Time_circadian_clock, ZT_initial)

Par_decay_cell = num2cell(Par_decay);
[d_CD4Tn, d_CD8Tn, d_Neu, d_APCi, d_APCm, d_Bn, d_BGC, d_Bm, d_ASC, d_Th, d_CD4Ta, d_CD4Tm, d_CTL, d_CD8Ta, d_CD8Tm, ...
    d_Treg, d_TNFa, d_IL10, d_CCL2, d_IL6, d_IFNg, d_IL2, d_IL4, d_CCL1921, d_CXCL5, d_CXCL8, d_IL1b, d_Treg_APCm, ...
    d_Treg_Neu, d_Treg_T, d_Ab, d_Ag, d_Ad, d_S1PR1, d_CA, d_D, d_IL10, d_IL102, d_IL6, d_N, d_TNF, d_DC, d_HEC] = deal(Par_decay_cell{:});

Par_produce_cell = num2cell(Par_produce);
[p0, a_0_APCi, a_Ag_APCm, a_0_Neu, a_0_CD4Tn, a_CD4Tn_CD4Ta, a_CD4Tm_CD4Ta, a_CD4Ta_Th1, a_CD4Ta_Th2, a_CD4Ta_Th17, a_CD4Ta_Tfh, ...
    a_CD4Ta_Treg, a_CD4Ta_CD4Tm, a_0_CD8Tn, a_CD8Tn_CD8Ta, a_CD8Tm_CD8Ta, a_CD8Ta_CTL, a_CD8Ta_CD8Tm, a_0_Bn, a_Bn_BGC, a_Bm_BGC, ...
    a_BGC_ASC, a_BGC_Bm, a_0_CCL2, a_APCm_CCL2, a_0_CCL1921, a_APCm_CCL1921, a_0_CXCL5, a_Neu_CXCL5, a_APCm_CXCL8, a_S1PR1, a_APCm_TNFa, ...
    a_Th1_TNFa, a_CTL_TNFa, a_Neu_IFNg, a_Th1_IFNg, a_CTL_IFNg, a_APCm_IL1b, a_Th1_IL2, a_CD4Ta_IL2, a_CD8Ta_IL2, a_Th2_IL4, a_APCm_IL6, ...
    a_Th17_IL6, a_APCm_IL10, a_Treg_IL10, a_ASC_Ab, a_Bm_Ab, a_CA, a_D, a_IL10, a_IL102, a_IL6, a_N, a_P, a_TNF, a_IL10D, a_DC, a_TNF_DC, a_CD4Ta_Th, ...
    ] = deal(Par_produce_cell{:});

Par_transition_cell = num2cell(Par_transition);
[c_hHEC_iHEC, c_iHEC_hHEC, c_PTB_cyto, c_BPT_cyto, c_PTLN_cyto, c_LNB_cyto, c_BLN_cyto, c_BPT_cell, c_PTLN_cell, c_LNB_cell, c_BLN_cell] = deal(Par_transition_cell{:});

Par_halfmax_cell = num2cell(Par_halfmax);
[K_TNFa_APCi, K_CCL2_APCi, K_IFNg_APCi, K_IL10_APCi, K_TNFa_APCm, K_CCL1921_APCm, K_TNFa_Neu, K_CXCL5_Neu, K_CXCL8_Neu, K_TNFa_CD4Tn, ...
    K_CCL1921_CD4Tn, K_S1PR1_CD4Tn, K_IL2_CD4Tn, K_APCm_CD4Tn, K_IL2_CD4Tm, K_APCm_CD4Tm, K_APCm_CD4Ta, K_IFNg_Th1, K_IL4_Th1, K_IL10_Th1, ...
    K_TNFa_Th1, K_CCL1921_Th1, K_S1PR1_Th1, K_IFNg_Th2, K_IL10_Th2, K_CCL1921_Th2, K_S1PR1_Th2, K_IL4_Th2, K_TNFa_Th2, K_IL1b_Th17, K_IL6_Th17, ...
    K_IL10_Th17, K_TNFa_Th17, K_CCL1921_Th17, K_S1PR1_Th17, K_IL10_Tfh, K_IL2_Treg, K_IL6_Treg, K_IL10_Treg, K_TNFa_Treg, K_CCL1921_Treg, ...
    K_S1PR1_Treg, K_TNFa_CD8Tn, K_CCL1921_CD8Tn, K_S1PR1_CD8Tn, K_IL2_CD8Tn, K_APCm_CD8Tn, K_IL2_CD8Tm, K_APCm_CD8Tm, K_APCm_CD8Ta, ...
    K_IL2_CTL, K_IL6_CTL, K_Th1_CTL, K_IL10_CTL, K_TNFa_CTL, K_CCL1921_CTL, K_S1PR1_CTL, K_TNFa_Bn, K_CCL1921_Bn, K_S1PR1_Bn, K_IL4_BGC, ...
    K_IL6_BGC, K_APCm_BGC, K_Tfh_BGC, K_TNFa_ASC, K_CCL1921_ASC, K_S1PR1_ASC, K_IL1b_CCL2, K_IL10_CCL2, K_IL1b_CXCL5, K_IL10_CXCL5, ...
    K_IL1b_CXCL8, K_IL10_CXCL8, K_IL1b_TNFa, K_IL10_TNFa, K_IL10_IL1b, K_IL1b_IL6, K_IL10_IL6, K_IL4_Ab, K_D, K_IL10, K_IL102, K_IL10d, ...
    K_IL10IL6, K_IL10REV, K_IL10TNF, K_IL6, K_IL6CA, K_IL6IL10, K_IL6IL6, K_IL6REV, K_IL6TNF, K_N, K_NCA, K_NIL10, K_NIL6, K_NTNF, K_P, ...
    K_TNFCA, K_TNFCRY, K_TNFIL10, K_TNFIL6, K_TNFROR, K_TNFTNF, K_CCL21, K_S1PR1, K_DC, K_DCCA, K_DCIL10, K_DCIFN, K_TNFa_Th, K_DCTNF, ...
    K_TNF_DC, K_IL2_CD4Ta, K_IL2_Th, K_IL6_Th, K_IFNg_Th, K_IL10_Th, K_APCm_iHEC, K_TNFa_iHEC, K_DC_CCL1921, K_IL2_CD8Ta, K_Th_CTL, K_DC_t, K_CCL19, ...
    K_CCL19_DC, K_CCL21_DC] = deal(Par_halfmax_cell{:});

Par_Others_cell = num2cell(Par_Others);
[N_GC, rho_BGC, V_B, V_PT, V_LN, s_CA, s_IL10, N_CD4Tn, N_CD8Tn, N_Bn, rho_CD4Ta, rho_CD8Ta, N_CD4Ta, N_CD8Ta, ...
    N_Neu, rho_HEC, N_HEC, sigma_CD4T, n_CCL1921_APCm, control_CRY_TNF, control_ROR_TNF, control_REV_IL6, control_REV_IL10, ...
    control_CB_CCL21, control_CB_S1PR1, control_PER, control_CRY, control_REV, control_ROR, control_BMAL1, ...
    control_DC_rhythmic_homing, control_mphi, control_EC, control_T, control_CpG, control_Lymphocyte_rhythmic_homing, ...
    control_cytokine_secretion, control_DCT_interaction, control_TB_proliferation, phi_deviation, control_bistable_innate, ...
    control_DC_migration_pre, lambda_crytnf, lambda_rortnf, lambda_revil6, lambda_revil10, lambda_cbccl21, lambda_cbs1pr1] = deal(Par_Others_cell{:});

Par_promotion_cell = num2cell(Par_promotion);
[alpha_TNFa_APCi, alpha_IFNg_APCi, alpha_TNFa_APCm, alpha_TNFa_Neu, alpha_TNFa_CD4Tn, alpha_IL2_CD4Tn, alpha_IL2_CD4Tm, alpha_IFNg_Th1, ...
    alpha_TNFa_Th1, alpha_IL4_Th2, alpha_TNFa_Th2, alpha_IL1b_Th17, alpha_IL6_Th17, alpha_TNFa_Th17, alpha_IL2_Treg, alpha_IL10_Treg, ...
    alpha_TNFa_Treg, alpha_TNFa_CD8Tn, alpha_IL2_CD8Tn, alpha_IL2_CD8Tm, alpha_IL2_CTL, alpha_IL6_CTL, alpha_Th1_CTL, alpha_TNFa_CTL, ...
    alpha_TNFa_Bn, alpha_IL4_BGC, alpha_IL6_BGC, alpha_Tfh_BGC, alpha_TNFa_ASC, alpha_IL1b_CCL2, alpha_IL1b_CXCL5, alpha_IL1b_CXCL8, ...
    alpha_IL1b_TNFa, alpha_IL1b_IL6, alpha_IL4_Ab, alpha_IL10IL6, alpha_IL10TNF, alpha_IL6IL6, alpha_IL6TNF, alpha_ND, alpha_NIL6, alpha_NP...
    alpha_NTNF, alpha_TNFTNF, alpha_DCP, alpha_DCTNF, alpha_DCIFN, alpha_TNFa_Th, alpha_TB_hHEC, alpha_TB_iHEC, alpha_IL2_CD4Ta, alpha_IL2_Th, ...
    alpha_IL6_Th, alpha_IFNg_Th, alpha_iHEC, alpha_APCm_iHEC, alpha_IL2_CD8Ta, alpha_Th_CTL, alpha_CCL19, alpha_CCL19_DC...
    ] = deal(Par_promotion_cell{:});

Par_Clock_immune_cell = num2cell(Par_Clock_immune);
[K_BMAL1_APCi, K_REV_CCL2, K_ROR_TNFa, K_REV_IL1b, K_REV_IL6, K_REV_IL10, K_REV_Th, K_BC_CCL1921, K_BC_S1PR1, ...
    K_IL1b_REV, alpha_IL1b_REV, n_IL1b_REV, K_BC_CCL2, K_BC_APC] = deal(Par_Clock_immune_cell{:});

gamma_BPT = V_B / V_PT; gamma_PTLN = V_PT / V_LN; gamma_LNB = V_LN / V_B; gamma_BLN = V_B / V_LN; gamma_PTB = V_PT / V_B;

yi_cell = num2cell(y(1 : 72));
[Ag_PT, Ad_PT, DC_PT, mphi_PT, Neu_PT, Treg_PT, Th_PT, Th2_PT, Th17_PT, CTL_PT, ASC_PT, CCL2_PT, ...   % 1~12
    CXCL5_PT, CXCL8_PT, TNFa_PT, IFNg_PT, IL1b_PT, IL2_PT, IL4_PT, IL6_PT, IL10_PT, Ab_PT, ...     % 13~22
    DC_LN, CD4Tn_LN, CD4Ta_LN, Th_LN, hHEC, iHEC, CCL1921_DC_LN, Treg_LN, CD4Tm_LN, ...           % 23~31
    CD8Tn_LN, CD8Ta_LN, CTL_LN, CD8Tm_LN, Bn_LN, BGC_LN, ASC_LN, Bm_LN, ...                         % 32~39
    siCCL21_LN, S1PR1, TNFa_LN, IFNg_LN, CA_LN, IL2_LN, IL4_LN, IL6_LN, IL10_LN, Ab_LN, ...       % 40~49
    sCCL21_LN, Neu_B, CD4Tn_B, CD8Tn_B, Bn_B, Th_B, Th2_B, Th17_B, Treg_B, CTL_B, ASC_B, ...          % 50~60
    TNFa_B, IFNg_B, IL1b_B, IL2_B, IL4_B, IL6_B, IL10_B, Ab_B, CA_B, D, CA_PT, Y_IL10] = deal(yi_cell{:}); % 61~68

% variables of clock model
yc_cell = num2cell(y(73 : 84));
[Per, Cry, Rev, Ror, Bmal1, PER, CRY, REV, ROR, ...
    BMAL1, PER_CRY, CLOCK_BMAL1] = deal(yc_cell{:});


% parameters
% mRNA and protein degradation rate constants (in h^-1)
Par_clock(1 : 12) = Par_clock(1 : 12) * 24;
par1 = num2cell(Par_clock(1 : 12));
[dm_per, dm_cry, dm_rev, dm_ror, dm_bmal, dp_per, ...
    dp_cry, dp_rev, dp_ror, dp_bmal, d_pc, d_cb] = deal(par1{:});

% maximal transcription rates (in nmol l^-1 h^-1)
Par_clock(13 : 17) = Par_clock(13 : 17) * 24;
par2 = num2cell(Par_clock(13 : 17));
[vmax_per, vmax_cry, vmax_rev, vmax_ror, vmax_bmal] = deal(par2{:});

% activation ratios (demensionless)
par3 = num2cell(Par_clock(18 : 22));
[fold_per, fold_cry, fold_rev, fold_ror, fold_bmal] = deal(par3{:});

% regulation thresholds (in nmol/l)
par4 = num2cell(Par_clock(23 : 33));
[Ka_per_cb, Ki_per_pc, Ka_cry_cb, Ki_cry_pc, Ki_cry_rev, Ka_rev_cb, ...
    Ki_rev_pc, Ka_ror_cb, Ki_ror_pc, Ka_bmal_ror, Ki_bmal_rev] = deal(par4{:});

% hill coefficients (dimensionless)
par5 = num2cell(Par_clock(34 : 44));
[hill_per_cb, hill_per_pc, hill_cry_cb, hill_cry_pc, ...
    hill_cry_rev, hill_rev_cb, hill_rev_pc, hill_ror_cb, ...
    hill_ror_pc, hill_bmal_ror, hill_bmal_rev] = deal(par5{:});

% translation rates (in molecules per hour per mRNA)
Par_clock(45 : 49) = Par_clock(45 : 49) * 24;
par6 = num2cell(Par_clock(45 : 49));
[kp_per, kp_cry, kp_rev, kp_ror, kp_bmal] = deal(par6{:});

% complexation kinetic rates
Par_clock(50 : 53) = Par_clock(50 : 53) * 24;
par7 = num2cell(Par_clock(50 : 53));
[kass_cb, kass_pc, kdiss_cb, kdiss_pc] = deal(par7{:});

% Purpose of introducing the circa_ periodic function: To achieve specific knockout
% During steady state, ZT_initial is uniformly set to 0
% In response phase, ZT_initial takes the ZT at the time of injection
% To ensure t ≥ 0, it is recommended to set ZT0 to 24 to avoid errors
% PER-1 CRY-2 REV_ERB-3 ROR-4 BMAL1_5 PER:CRY-6 CLOCK:BMAL1-7
circa_t_index = mod(floor((t * 24 + ZT_initial) / 24 * 1000), 1000) + 1;
circa_t_dev_index = mod(floor((t * 24 + ZT_initial + phi_deviation) / 24 * 1000), 1000) + 1;
% PER_t=Time_circadian_clock(circa_t_index,1);
CRY_t = Time_circadian_clock(circa_t_index, 2);
REV_t = Time_circadian_clock(circa_t_index, 3);
ROR_t = Time_circadian_clock(circa_t_index, 4);
% BMAL1_t=Time_circadian_clock(circa_t_index,5);
% PER_CRY_t=Time_circadian_clock(circa_t_index,6);
CLOCK_BMAL1_t = Time_circadian_clock(circa_t_dev_index, 7);
% par_cyto_log_all = zeros(1, 53);
% par_cyto_log_all(2 : 52) = par_cyto_log_1;
% par_cyto_1 = par_0_1 .* 10 .^ par_cyto_log_all;
% par_cyto_1_cell = num2cell(par_cyto_1);
% [d_Ag, a_N, K_N, d_N, alpha_NP, alpha_ND, K_NTNF, K_NIL6, K_NCA, K_NIL10, alpha_NTNF, alpha_NIL6, ...
%     a_D, d_D, K_D, a_P, K_P, ...
%     a_CA, d_CA, s_CA, ...
%     alpha_IL6TNF, K_IL6TNF, a_IL6, d_IL6, K_IL6, K_IL6IL10, alpha_IL6IL6, K_IL6IL6, K_IL6CA, ...
%     a_TNF, d_TNF, K_TNFIL10, K_TNFCA, alpha_TNFTNF, K_TNFTNF, K_TNFIL6, ...
%     alpha_IL10TNF, K_IL10TNF, alpha_IL10IL6, K_IL10IL6, a_IL10, d_IL10, K_IL10, s_IL10, K_IL10d, ...
%     a_IL102, d_IL102, K_IL102, K_IL6REV, K_IL10REV, K_TNFROR, K_TNFCRY, p0] = par_cyto_1_cell{:};

n_S1PR1 = 3;
n_sec_CCL19 = 2;
n_CCL1921_APCm = 2;
n_DC_HEC = 1;
K_DC = 1;
a = 1;
b = 1;
a_DC = 205;
K_IL10d = 5.3 * 10 ^ 12;
K_YIL10 = 7 * 10 ^ 5;
% rho_BGC = 4;
% a_BGC_ASC=100*a_BGC_ASC;
% a_ASC_Ab = 0.1*a_ASC_Ab;
K = 0.15;
% K_DC = 0.01;
% a_DC = 15;


% a_IL10D = 0;
if control_CpG == 1
    a_TNF = 0.1 * a_TNF;
    %     a_IL10 = 0;
    a_IL6 = 0.1 * a_IL6;
    K_DC = 4000;
    a_DC = 200 * a_DC;
    %     alpha_DCP=alpha_DCP;
    d_Ag = 72;
    %     a_IL102 = 0;
end

dydt = zeros(84, 1);
% Antigen
dydt(1) = - d_Ag * Ag_PT;
% This variable is no longer used
dydt(2) = - d_Ad * Ad_PT;

%DC
if control_CpG == 1
    R_DC = alpha_DCP * Ag_PT * 1 / (1 + (CA_PT / K_DCCA)) * 1 / (1 + (IL10_PT / K_DCIL10)) * 2 / (1 + (CLOCK_BMAL1_t / K_BC_APC) ^ 3) * (1 + alpha_DCTNF * TNFa_PT / (TNFa_PT + K_DCTNF )...
        + alpha_DCIFN * IFNg_PT / (IFNg_PT + K_DCIFN ));
    dydt(3) = a_DC * R_DC / (K_DC + R_DC) - c_PTLN_cell * b * (CCL1921_DC_LN ^ n_CCL1921_APCm) / (K_CCL1921_APCm ^ n_CCL1921_APCm + (CCL1921_DC_LN) ^ n_CCL1921_APCm) * DC_PT  ...
        - d_DC * DC_PT;
end
if control_CpG == 0

    R_DC = alpha_DCP * Ag_PT * 1 / (1 + (CA_PT / K_DCCA)) * 1 / (1 + (IL10_PT / K_DCIL10)) * (1 + alpha_DCTNF * TNFa_PT / (TNFa_PT + K_DCTNF )...
        + alpha_DCIFN * IFNg_PT / (IFNg_PT + K_DCIFN ));

    dydt(3) = a_DC * R_DC / (K_DC + R_DC) - c_PTLN_cell * b * (CCL1921_DC_LN ^ n_CCL1921_APCm) / (K_CCL1921_APCm ^ n_CCL1921_APCm + (CCL1921_DC_LN) ^ n_CCL1921_APCm) * DC_PT  ...
        - d_DC * DC_PT;

end


%mphi
if control_CpG == 0
    R_mphi = (alpha_NP * Ag_PT + alpha_ND * D) * 1 / (1 + (CA_PT / K_NCA)) * 1 / (1 + (IL10_PT / K_NIL10)) * (1 + alpha_NTNF * TNFa_PT / (TNFa_PT + K_NTNF )) * (1 + alpha_NIL6 * IL6_PT / (IL6_PT + K_NIL6 ));
    dydt(4) = a_N * R_mphi / (K_N + R_mphi) - d_N * mphi_PT;
end
if control_CpG == 1
    R_mphi = (alpha_NP * Ag_PT + alpha_ND * D) * 1 / (1 + (CA_PT / K_NCA)) * 1 / (1 + (IL10_PT / K_NIL10)) * 2 / (1 + (CLOCK_BMAL1_t / K_BC_APC) ^ 3) * (1 + alpha_NTNF * TNFa_PT / (TNFa_PT + K_NTNF )) * (1 + alpha_NIL6 * IL6_PT / (IL6_PT + K_NIL6 ));
    dydt(4) = a_N * R_mphi / (K_N + R_mphi) - d_N * mphi_PT;
end

%Neu
dydt(5) = gamma_BPT * c_BPT_cell * (1 + alpha_TNFa_Neu * TNFa_PT / (K_TNFa_Neu + TNFa_PT)) * (CXCL5_PT ^ 2 / (K_CXCL5_Neu ^ 2 + CXCL5_PT ^ 2) + CXCL8_PT ^ 2 / (K_CXCL8_Neu ^ 2 + CXCL8_PT ^ 2)) * Neu_B - ...
    d_Neu * Neu_PT - d_Treg_Neu * Neu_PT * Treg_PT;

%Treg
dydt(6) = gamma_BPT * c_BPT_cell * (1 + alpha_TNFa_Treg * TNFa_PT / (K_TNFa_Treg + TNFa_PT)) * Treg_B - c_PTLN_cell * Treg_PT - d_Treg * Treg_PT;

%Th
dydt(7) = gamma_BPT * c_BPT_cell * (1 + alpha_TNFa_Th * TNFa_PT / (K_TNFa_Th + TNFa_PT)) * Th_B - c_PTLN_cell * Th_PT - d_Th * Th_PT - d_Treg_T * Th_PT * Treg_PT;

% dydt(7) = gamma_BPT * c_BPT_Th1 * (1 + alpha_TNFa_Th1 * TNFa_PT / (K_TNFa_Th1 + TNFa_PT)) * Th1_B - c_PTLN_Th1 * Th1_PT - d_Th * Th1_PT - d_Treg_T * Th1_PT * Treg_PT;

% This variable is no longer used
dydt(8) = 0;

% This variable is no longer used
dydt(9) = 0;

% CTL
dydt(10) = gamma_BPT * c_BPT_cell * (1 + alpha_TNFa_CTL * TNFa_PT / (K_TNFa_CTL + TNFa_PT)) * CTL_B - c_PTLN_cell * CTL_PT - d_CTL * CTL_PT - d_Treg_T * CTL_PT * Treg_PT;

% ASC
dydt(11) = gamma_BPT * c_BPT_cell * (1 + alpha_TNFa_ASC * TNFa_PT / (K_TNFa_ASC + TNFa_PT)) * ASC_B - c_PTLN_cell * ASC_PT - d_ASC * ASC_PT;

% CCL2(This variable is obsolete)
dydt(12) = 0;

% CXCL5
dydt(13) = a_0_CXCL5 + a_Neu_CXCL5 * K_IL10_CXCL5 / ( IL10_PT + K_IL10_CXCL5) * Neu_PT - d_CXCL5 * CXCL5_PT;
%  dydt(13) = 0;
% CXCL8
dydt(14) = a_APCm_CXCL8 * K_IL10_CXCL8 / ( IL10_PT + K_IL10_CXCL8) * mphi_PT - d_CXCL8 * CXCL8_PT;


% TNFa
if control_mphi == 1     
    if control_cytokine_secretion == 1
        dydt(15) = a_TNF * mphi_PT ^ 1.5 * (1 + alpha_TNFTNF * TNFa_PT / (TNFa_PT + K_TNFTNF )) * 1 / (1 + (CA_PT / K_TNFCA) ^ 6) * 1 / (1 + (IL10_PT / K_TNFIL10)) * 1 / (1 + (IL6_PT / K_TNFIL6)) * ...
            (1 / (1 + (CRY_t / K_TNFCRY))) * (1 / (1 + (ROR_t / K_TNFROR))) + a_TNF_DC * DC_PT ^ 4 / (DC_PT ^ 4 + K_TNF_DC ^ 4) + ...
            a_Th1_TNFa * Th_PT + a_CTL_TNFa * CTL_PT + gamma_BPT * c_BPT_cyto * TNFa_B - c_PTLN_cyto * TNFa_PT - c_PTB_cyto * TNFa_PT - d_TNF * TNFa_PT;
    else
        dydt(15) = a_TNF * mphi_PT ^ 1.5 * (1 + alpha_TNFTNF * TNFa_PT / (TNFa_PT + K_TNFTNF )) * 1 / (1 + (CA_PT / K_TNFCA) ^ 6) * 1 / (1 + (IL10_PT / K_TNFIL10)) * 1 / (1 + (IL6_PT / K_TNFIL6)) * ...
            0.471248991032103 * 0.404997065539742 + a_TNF_DC * DC_PT ^ 4 / (DC_PT ^ 4 + K_TNF_DC ^ 4) + ...
            a_Th1_TNFa * Th_PT + a_CTL_TNFa * CTL_PT + gamma_BPT * c_BPT_cyto * TNFa_B - c_PTLN_cyto * TNFa_PT - c_PTB_cyto * TNFa_PT - d_TNF * TNFa_PT;
    end
else
    dydt(15) = a_TNF * mphi_PT ^ 1.5 * (1 + alpha_TNFTNF * TNFa_PT / (TNFa_PT + K_TNFTNF )) * 1 / (1 + (CA_PT / K_TNFCA) ^ 6) * 1 / (1 + (IL10_PT / K_TNFIL10)) * 1 / (1 + (IL6_PT / K_TNFIL6)) * ...
        (0.461721323663885 + lambda_crytnf * (1 / (1 + (CRY_t / K_TNFCRY)) - 0.461721323663885)) * ( 0.404802698918162 + lambda_rortnf * (1 / (1 + (ROR / K_TNFROR)) - 0.404802698918162 )) + a_TNF_DC * DC_PT ^ 4 / (DC_PT ^ 4 + K_TNF_DC ^ 4) + ...
        a_Th1_TNFa * Th_PT + a_CTL_TNFa * CTL_PT + gamma_BPT * c_BPT_cyto * TNFa_B - c_PTLN_cyto * TNFa_PT - c_PTB_cyto * TNFa_PT - d_TNF * TNFa_PT; % 否则mphi使用动态行为下的生物钟（敲除，0态出发）
end

% IFNg
dydt(16) = a_Neu_IFNg * Neu_PT + a_CTL_IFNg * CTL_PT + gamma_BPT * c_BPT_cyto * IFNg_B - c_PTLN_cyto * IFNg_PT - c_PTB_cyto * IFNg_PT - d_IFNg * IFNg_PT;

% IL1b.This variable is no longer used
dydt(17) = 0;

% IL2
dydt(18) = a_Th1_IL2 * Th_PT + gamma_BPT * c_BPT_cyto * IL2_B - c_PTLN_cyto * IL2_PT - c_PTB_cyto * IL2_PT - d_IL2 * IL2_PT;

% IL4
dydt(19) = a_Th2_IL4 * Th_PT + gamma_BPT * c_BPT_cyto * IL4_B - c_PTLN_cyto * IL4_PT - c_PTB_cyto * IL4_PT - d_IL4 * IL4_PT;

% IL6
if control_mphi == 1     
    if control_cytokine_secretion == 1
        dydt(20) = a_IL6 * mphi_PT ^ 4 / (mphi_PT ^ 4 + K_IL6 ^ 4) * 1 / (1 + (IL10_PT / K_IL6IL10)) * 1 / (1 + (CA_PT / K_IL6CA)) * (1 / (1 + (REV_t / K_IL6REV)) ) * ...
            (1 + alpha_IL6TNF * TNFa_PT / (TNFa_PT + K_IL6TNF) + alpha_IL6IL6 * IL6_PT / (IL6_PT + K_IL6IL6)) + a_Th17_IL6 * Th_PT + gamma_BPT * c_BPT_cyto * IL6_B - ...
            c_PTLN_cyto * IL6_PT - c_PTB_cyto * IL6_PT - d_IL6 * IL6_PT;
    else
        dydt(20) = a_IL6 * mphi_PT ^ 4 / (mphi_PT ^ 4 + K_IL6 ^ 4) * 1 / (1 + (IL10_PT / K_IL6IL10)) * 1 / (1 + (CA_PT / K_IL6CA)) * 0.499069040924711 * ...
            (1 + alpha_IL6TNF * TNFa_PT / (TNFa_PT + K_IL6TNF) + alpha_IL6IL6 * IL6_PT / (IL6_PT + K_IL6IL6)) + a_Th17_IL6 * Th_PT + gamma_BPT * c_BPT_cyto * IL6_B - ...
            c_PTLN_cyto * IL6_PT - c_PTB_cyto * IL6_PT - d_IL6 * IL6_PT;
    end
else
    dydt(20) = a_IL6 * mphi_PT ^ 4 / (mphi_PT ^ 4 + K_IL6 ^ 4) * 1 / (1 + (IL10_PT / K_IL6IL10)) * 1 / (1 + (CA_PT / K_IL6CA)) * (0.457178730750814 + lambda_revil6 * (1 / (1 + (REV / K_IL6REV)) - 0.457178730750814)) * ...
        (1 + alpha_IL6TNF * TNFa_PT / (TNFa_PT + K_IL6TNF) + alpha_IL6IL6 * IL6_PT / (IL6_PT + K_IL6IL6)) + a_Th17_IL6 * Th_PT + gamma_BPT * c_BPT_cyto * IL6_B - c_PTLN_cyto * IL6_PT - c_PTB_cyto * IL6_PT - d_IL6 * IL6_PT; % 否则mphi使用动态行为下的生物钟（敲除，0态出发）
end


% IL10
if control_mphi == 1     
    if control_cytokine_secretion == 1
        dydt(21) = a_IL10 * mphi_PT ^ 3 / (mphi_PT ^ 3 + K_IL10 ^ 3) * (1 + alpha_IL10IL6 * IL6_PT ^ 4 / (IL6_PT ^ 4 + K_IL10IL6 ^ 4)...
            + alpha_IL10TNF * TNFa_PT / (TNFa_PT + K_IL10TNF )) * (1 / (1 + (REV_t / K_IL10REV) )) + a_Treg_IL10 * Treg_PT + ...
            gamma_BPT * c_BPT_cyto * IL10_B - c_PTLN_cyto * IL10_PT - c_PTB_cyto * IL10_PT - d_IL10 / (1 + (IL10_PT / K_IL10d)) * IL10_PT + a_IL10D * (K_YIL10 / (K_YIL10 + Y_IL10)) * Y_IL10 + s_IL10;
    else
        dydt(21) = (a_IL10 * mphi_PT ^ 3 / (mphi_PT ^ 3 + K_IL10 ^ 3) * (1 + alpha_IL10IL6 * IL6_PT ^ 4 / (IL6_PT ^ 4 + K_IL10IL6 ^ 4)...
            + alpha_IL10TNF * TNFa_PT / (TNFa_PT + K_IL10TNF )) * 0.413180898459689 + a_Treg_IL10 * Treg_PT + ...
            gamma_BPT * c_BPT_cyto * IL10_B - c_PTLN_cyto * IL10_PT - c_PTB_cyto * IL10_PT - d_IL10 / (1 + (IL10_PT / K_IL10d)) * IL10_PT + a_IL10D * (K_YIL10 / (K_YIL10 + Y_IL10)) * Y_IL10 + s_IL10);
    end
else
    dydt(21) = a_IL10 * mphi_PT ^ 3 / (mphi_PT ^ 3 + K_IL10 ^ 3) * (1 + alpha_IL10IL6 * IL6_PT ^ 4 / (IL6_PT ^ 4 + K_IL10IL6 ^ 4)...
        + alpha_IL10TNF * TNFa_PT / (TNFa_PT + K_IL10TNF )) * (0.361645774990316 + lambda_revil10 * (1 / (1 + (REV / K_IL10REV)) - 0.361645774990316 )) + a_Treg_IL10 * Treg_PT + ...
        gamma_BPT * c_BPT_cyto * IL10_B - c_PTLN_cyto * IL10_PT - c_PTB_cyto * IL10_PT - d_IL10 / (1 + (IL10_PT / K_IL10d)) * IL10_PT + a_IL10D * (K_YIL10 / (K_YIL10 + Y_IL10)) * Y_IL10 + s_IL10; % 否则mphi使用动态行为下的生物钟（敲除，0态出发）
end

% Antibody
dydt(22) = a_ASC_Ab * (1 + alpha_IL4_Ab * IL4_PT / (IL4_PT + K_IL4_Ab)) * ASC_PT + gamma_BPT * c_BPT_cyto * Ab_B - c_PTLN_cyto * Ab_PT - c_PTB_cyto * Ab_PT - d_Ab * Ab_PT;

alpha_HEC = alpha_TB_hHEC * hHEC + alpha_TB_iHEC * iHEC;
alpha_HEC_CD4T = alpha_TB_hHEC * hHEC + sigma_CD4T * alpha_TB_iHEC * iHEC;
alpha_CCL1921 = siCCL21_LN ^ 2 / (K_CCL21 ^ 2 + siCCL21_LN ^ 2) + alpha_CCL19 * sCCL21_LN ^ 2 / (K_CCL19 ^ 2 + sCCL21_LN ^ 2);
alpha_CCL1921_DC = siCCL21_LN ^ 2 / (K_CCL21_DC ^ 2 + siCCL21_LN ^ 2) + alpha_CCL19_DC * sCCL21_LN ^ 2 / (K_CCL19_DC ^ 2 + sCCL21_LN ^ 2);
alpha_egress = S1PR1 ^ n_S1PR1 / (S1PR1 ^ n_S1PR1 + K_S1PR1 ^ n_S1PR1);

% DC_LN

dydt(23) = gamma_PTLN * c_PTLN_cell * alpha_CCL1921_DC * DC_PT - d_DC * DC_LN;


if control_DCT_interaction == 1
    alpha_activation_CB = CLOCK_BMAL1_t / (CLOCK_BMAL1_t + 0.1649);
else
    alpha_activation_CB = 0.486344862333709;
end
DI7 = 0.93; Lambda = 0.27; phi_p = 4 / 24; omega = 2 * pi;
if control_TB_proliferation == 1
    DI = DI7 + Lambda * (1 + sin(omega * (t + ZT_initial / 24 + phi_p + phi_deviation / 24)));
else
    DI = 1.2;
end
% CD4Tn
dydt(24) = gamma_BLN * c_BLN_cell * alpha_HEC_CD4T * alpha_CCL1921 * CD4Tn_B...
    - c_LNB_cell * alpha_egress * CD4Tn_LN - a_CD4Tn_CD4Ta * alpha_activation_CB * (1 + alpha_IL2_CD4Tn * IL2_LN / (K_IL2_CD4Tn + IL2_LN)) * ...
    DC_LN / (DC_LN + K_APCm_CD4Tn) * CD4Tn_LN - d_CD4Tn * CD4Tn_LN;

% CD4Ta
dydt(25) = a_CD4Tn_CD4Ta * alpha_activation_CB * (1 + alpha_IL2_CD4Tn * IL2_LN / (K_IL2_CD4Tn + IL2_LN)) * DC_LN / (DC_LN + K_APCm_CD4Tn) * CD4Tn_LN + ...
    a_CD4Tm_CD4Ta * (1 + alpha_IL2_CD4Tm * IL2_LN / (K_IL2_CD4Tm + IL2_LN)) * DC_LN / (DC_LN + K_APCm_CD4Tm) * CD4Tm_LN + ...
    rho_CD4Ta * DI * CD4Ta_LN * (1 - CD4Ta_LN / N_CD4Ta ) * (1 + alpha_IL2_CD4Ta * IL2_LN / (K_IL2_CD4Ta + IL2_LN)) * DC_LN / (DC_LN + K_APCm_CD4Tn)...
    - d_CD4Ta * CD4Ta_LN - d_Treg_T * CD4Ta_LN * Treg_LN;


% Th

dydt(26) = a_CD4Ta_Th * 2 * REV / (REV + K_REV_Th) * (1 + alpha_IL2_Th * IL2_LN / (IL2_LN + K_IL2_Th) + alpha_IL6_Th * IL6_LN / (IL6_LN + K_IL6_Th) + alpha_IFNg_Th * IFNg_LN / (IFNg_LN + K_IFNg_Th)) * ...
    K_IL10_Th / (K_IL10_Th + IL10_LN) * K ^ 2 / (K ^ 2 + DC_LN ^ 2) * CD4Ta_LN - d_Treg_T * Th_LN * Treg_LN + gamma_BLN * alpha_HEC_CD4T * c_BLN_cell * alpha_CCL1921 * Th_B...
    - c_LNB_cell * alpha_egress * Th_LN - d_Th * Th_LN + gamma_PTLN * c_PTLN_cell * Th_PT;



% hHEC

dydt(27) = -c_hHEC_iHEC * (DC_LN ^ n_DC_HEC / (DC_LN ^ n_DC_HEC + K_APCm_iHEC ^ n_DC_HEC) + alpha_iHEC * TNFa_LN / (TNFa_LN + K_TNFa_iHEC)) * hHEC + ...
    c_iHEC_hHEC * iHEC + rho_HEC * hHEC * (1 - (hHEC + iHEC) / N_HEC) ...
    - d_HEC * hHEC;
%iHEC

dydt(28) = c_hHEC_iHEC * (DC_LN ^ n_DC_HEC / (DC_LN ^ n_DC_HEC + K_APCm_iHEC ^ n_DC_HEC) + alpha_iHEC * TNFa_LN / (TNFa_LN + K_TNFa_iHEC)) * hHEC + ...
    rho_HEC * iHEC * (1 - (hHEC + iHEC) / N_HEC) ...
    - c_iHEC_hHEC * iHEC - d_HEC * iHEC;

% CCL1921_DC. This variable is no longer used
if control_DC_rhythmic_homing == 1
    dydt(29) = a_0_CCL1921 * (0.473053544293661 + lambda_cbccl21 * (CLOCK_BMAL1 .^ 3 ./ (K_BC_CCL1921 .^ 3 + CLOCK_BMAL1 .^ 3) - 0.473053544293661)) + ...
        a * a_APCm_CCL1921 * DC_LN ^ n_sec_CCL19 / (DC_LN ^ n_sec_CCL19 + K_DC_CCL1921 ^ n_sec_CCL19) - d_CCL1921 * CCL1921_DC_LN;
else
    dydt(29) = a_0_CCL1921 * 0.467 + ...
        a * a_APCm_CCL1921 * DC_LN ^ n_sec_CCL19 / (DC_LN ^ n_sec_CCL19 + K_DC_CCL1921 ^ n_sec_CCL19) - d_CCL1921 * CCL1921_DC_LN;
end

% Treg
dydt(30) = a_CD4Ta_Treg * (1 + alpha_IL10_Treg * IL10_LN / (IL10_LN + K_IL10_Treg) + alpha_IL2_Treg * IL2_LN / (IL2_LN + K_IL2_Treg)) * ...
    K_IL6_Treg / (K_IL6_Treg + IL6_LN) * CD4Ta_LN + ...
    gamma_BLN * c_BLN_cell * alpha_HEC_CD4T * alpha_CCL1921 * Treg_B...
    - c_LNB_cell * alpha_egress * Treg_LN - d_Treg * Treg_LN + gamma_PTLN * c_PTLN_cell * Treg_PT;

% CD4Tm
dydt(31) = a_CD4Ta_CD4Tm * DC_LN / (K_APCm_CD4Ta + DC_LN) * CD4Ta_LN - ...
    a_CD4Tm_CD4Ta * (1 + alpha_IL2_CD4Tm * IL2_LN / (K_IL2_CD4Tm + IL2_LN)) * DC_LN / (DC_LN + K_APCm_CD4Tm) * CD4Tm_LN - ...
    d_CD4Tm * CD4Tm_LN;

% CD8Tn
dydt(32) = gamma_BLN * c_BLN_cell * alpha_HEC...
    * alpha_CCL1921 * CD8Tn_B...
    - c_LNB_cell * alpha_egress * CD8Tn_LN - a_CD8Tn_CD8Ta * alpha_activation_CB * (1 + alpha_IL2_CD8Tn * IL2_LN / (K_IL2_CD8Tn + IL2_LN)) * ...
    DC_LN / (DC_LN + K_APCm_CD8Tn) * CD8Tn_LN - d_CD8Tn * CD8Tn_LN;

% CD8Ta
dydt(33) = a_CD8Tn_CD8Ta * alpha_activation_CB * (1 + alpha_IL2_CD8Tn * IL2_LN / (K_IL2_CD8Tn + IL2_LN)) * DC_LN / (DC_LN + K_APCm_CD8Tn) * CD8Tn_LN + ...
    a_CD8Tm_CD8Ta * (1 + alpha_IL2_CD8Tm * IL2_LN / (K_IL2_CD8Tm + IL2_LN)) * DC_LN / (DC_LN + K_APCm_CD8Tm) * CD8Tm_LN + ...
    rho_CD8Ta * DI * CD8Ta_LN * (1 - CD8Ta_LN / N_CD8Ta ) * (1 + alpha_IL2_CD8Ta * IL2_LN / (K_IL2_CD8Ta + IL2_LN)) * DC_LN / (DC_LN + K_APCm_CD8Tn)...
    - d_CD8Ta * CD8Ta_LN - d_Treg_T * CD8Ta_LN * Treg_LN;

% CTL
dydt(34) = a_CD8Ta_CTL * (1 + alpha_IL2_CTL * IL2_LN / (K_IL2_CTL + IL2_LN) + alpha_IL6_CTL * IL6_LN / (K_IL6_CTL + IL6_LN) + alpha_Th_CTL * Th_LN / (K_Th_CTL + Th_LN))...
    * K_IL10_CTL / (K_IL10_CTL + IL10_LN) * K ^ 2 / (K ^ 2 + DC_LN ^ 2) * CD8Ta_LN + ...
    gamma_BLN * c_BLN_cell * alpha_HEC * alpha_CCL1921 * CTL_B...
    - c_LNB_cell * alpha_egress * CTL_LN - d_CTL * CTL_LN - d_Treg_T * CTL_LN * Treg_LN + gamma_PTLN * c_PTLN_cell * CTL_PT;

% CD8Tm
dydt(35) = a_CD8Ta_CD8Tm * DC_LN / (K_APCm_CD8Ta + DC_LN) * CD8Ta_LN - ...
    a_CD8Tm_CD8Ta * (1 + alpha_IL2_CD8Tm * IL2_LN / (K_IL2_CD8Tm + IL2_LN)) * DC_LN / (DC_LN + K_APCm_CD8Tm) * CD8Tm_LN - ...
    d_CD8Tm * CD8Tm_LN;

% Bn
dydt(36) = gamma_BLN * c_BLN_cell * alpha_HEC * alpha_CCL1921 * Bn_B...
    - c_LNB_cell * alpha_egress * Bn_LN - a_Bn_BGC * alpha_activation_CB * (1 + alpha_IL4_BGC * IL4_LN / (K_IL4_BGC + IL4_LN) + alpha_IL6_BGC * IL6_LN / (K_IL6_BGC + IL6_LN)) * ...
    DC_LN / (DC_LN + K_APCm_BGC) * Bn_LN - d_Bn * Bn_LN;

% BGC
dydt(37) = DC_LN / (DC_LN + K_APCm_BGC) * (1 + alpha_IL4_BGC * IL4_LN / (K_IL4_BGC + IL4_LN) + alpha_IL6_BGC * IL6_LN / (K_IL6_BGC + IL6_LN)) * ...
    (a_Bn_BGC * alpha_activation_CB * Bn_LN + a_Bm_BGC * Bm_LN + rho_BGC * DI * alpha_Tfh_BGC * Th_LN / (K_Tfh_BGC + Th_LN) * BGC_LN * (1 - BGC_LN / N_GC)) ...
    - d_BGC * BGC_LN;

% ASC
dydt(38) = gamma_BLN * c_BLN_cell * alpha_HEC * alpha_CCL1921 * ASC_B...
    - c_LNB_cell * alpha_egress * ASC_LN + gamma_PTLN * c_PTLN_cell * ASC_PT + a_BGC_ASC * BGC_LN * K ^ 2 / (K ^ 2 + DC_LN ^ 2)...
    - d_ASC * ASC_LN;

% Bm
dydt(39) = a_BGC_Bm * BGC_LN - a_Bm_BGC * (1 + alpha_IL4_BGC * IL4_LN / (K_IL4_BGC + IL4_LN) + alpha_IL6_BGC * IL6_LN / (K_IL6_BGC + IL6_LN)) * ...
    DC_LN / (DC_LN + K_APCm_BGC) * Bm_LN - d_Bm * Bm_LN;

% siCCL21
if control_Lymphocyte_rhythmic_homing == 1     
    dydt(40) = a_0_CCL1921 * (CLOCK_BMAL1 .^ 3 ./ (K_BC_CCL1921 .^ 3 + CLOCK_BMAL1 .^ 3) ) - d_CCL1921 * siCCL21_LN;
else
    dydt(40) = a_0_CCL1921 * 0.467 - d_CCL1921 * siCCL21_LN; 
end

% S1PR1
if control_Lymphocyte_rhythmic_homing == 1     
    dydt(41) = a_S1PR1 * (K_BC_S1PR1 .^ 3 ./ (K_BC_S1PR1 .^ 3 + (CLOCK_BMAL1_t .^ 3) )) - d_S1PR1 * S1PR1;
else
    dydt(41) = a_S1PR1 * 0.0634 - d_S1PR1 * S1PR1; 
end


% TNFa
dydt(42) = a_TNF_DC * DC_LN ^ 4 / (DC_LN ^ 4 + K_TNF_DC ^ 4) + a_Th1_TNFa * Th_LN + a_CTL_TNFa * CTL_LN + gamma_PTLN * c_PTLN_cyto * TNFa_PT + ...
    gamma_BLN * c_BLN_cyto * TNFa_B - c_LNB_cyto * TNFa_LN - d_TNF * TNFa_LN;

% IFNg
dydt(43) = a_CTL_IFNg * CTL_LN + gamma_PTLN * c_PTLN_cyto * IFNg_PT + gamma_BLN * c_BLN_cyto * IFNg_B - c_LNB_cyto * IFNg_LN - d_IFNg * IFNg_LN;

% CA_LN
dydt(44) = - d_CA * CA_LN + gamma_PTLN * c_PTLN_cyto * CA_PT + gamma_BLN * c_BLN_cyto * CA_B - c_LNB_cyto * CA_LN;

% IL2
dydt(45) = a_Th1_IL2 * Th_LN + a_CD4Ta_IL2 * CD4Ta_LN + a_CD8Ta_IL2 * CD8Ta_LN + gamma_PTLN * c_PTLN_cyto * IL2_PT + ...
    gamma_BLN * c_BLN_cyto * IL2_B - c_LNB_cyto * IL2_LN - d_IL2 * IL2_LN;

% IL4
dydt(46) = a_Th2_IL4 * Th_LN + gamma_PTLN * c_PTLN_cyto * IL4_PT + gamma_BLN * c_BLN_cyto * IL4_B - c_LNB_cyto * IL4_LN - d_IL4 * IL4_LN;

% IL6
dydt(47) = a_Th17_IL6 * Th_LN + gamma_PTLN * c_PTLN_cyto * IL6_PT + gamma_BLN * c_BLN_cyto * IL6_B - c_LNB_cyto * IL6_LN - d_IL6 * IL6_LN;

% IL10
dydt(48) = (a_Treg_IL10 * Treg_LN + gamma_PTLN * c_PTLN_cyto * IL10_PT + gamma_BLN * c_BLN_cyto * IL10_B - c_LNB_cyto * IL10_LN - d_IL10...
    / (1 + (IL10_LN / K_IL10d)) * IL10_LN);

% Ab
dydt(49) = (a_ASC_Ab * ASC_LN + a_Bm_Ab * Bm_LN) * (1 + alpha_IL4_Ab * IL4_LN / (IL4_LN + K_IL4_Ab)) + gamma_PTLN * c_PTLN_cyto * Ab_PT + ...
    gamma_BLN * c_BLN_cyto * Ab_B - c_LNB_cyto * Ab_LN - d_Ab * Ab_LN;
% sCCL21
dydt(50) = a * a_APCm_CCL1921 * DC_LN ^ n_sec_CCL19 / (DC_LN ^ n_sec_CCL19 + K_DC_CCL1921 ^ n_sec_CCL19) - d_CCL1921 * sCCL21_LN;

% Neu(in Blood)
dydt(51) = a_0_Neu * Neu_B * (1 - Neu_B / N_Neu) - c_BPT_cell * (1 + alpha_TNFa_Neu * TNFa_PT / (K_TNFa_Neu + TNFa_PT))...
    * (CXCL5_PT ^ 2 / (K_CXCL5_Neu ^ 2 + CXCL5_PT ^ 2) + CXCL8_PT ^ 2 / (K_CXCL8_Neu ^ 2 + CXCL8_PT ^ 2)) * Neu_B - ...
    d_Neu * Neu_B;
% dydt(51) = 0;
% CD4Tn
dydt(52) = a_0_CD4Tn * CD4Tn_B * (1 - CD4Tn_B / N_CD4Tn) + gamma_LNB * c_LNB_cell * alpha_egress * CD4Tn_LN - ...
    c_BLN_cell * alpha_HEC_CD4T * alpha_CCL1921 * CD4Tn_B - d_CD4Tn * CD4Tn_B;

% CD8Tn
dydt(53) = a_0_CD8Tn * CD8Tn_B * (1 - CD8Tn_B / N_CD8Tn) + gamma_LNB * c_LNB_cell * alpha_egress * CD8Tn_LN - ...
    c_BLN_cell * alpha_HEC * alpha_CCL1921 * CD8Tn_B - d_CD8Tn * CD8Tn_B;

% Bn
dydt(54) = a_0_Bn * Bn_B * (1 - Bn_B / N_Bn) + gamma_LNB * c_LNB_cell * alpha_egress * Bn_LN - ...
    c_BLN_cell * alpha_HEC * alpha_CCL1921 * Bn_B - d_Bn * Bn_B;

% Th
dydt(55) = -c_BLN_cell * alpha_HEC_CD4T * alpha_CCL1921 * Th_B...
    + gamma_LNB * c_LNB_cell * alpha_egress * Th_LN - d_Th * Th_B - c_BPT_cell * (1 + alpha_TNFa_Th * TNFa_PT / (K_TNFa_Th + TNFa_PT)) * Th_B;

if control_DC_migration_pre == 1
    dydt(56) = gamma_PTLN * b * c_PTLN_cell * CCL1921_DC_LN ^ n_CCL1921_APCm / (K_CCL1921_APCm ^ n_CCL1921_APCm + CCL1921_DC_LN ^ n_CCL1921_APCm) * DC_PT;
else
    dydt(56) = 0;
end
dydt(57) = 0;

% Treg
dydt(58) = -c_BLN_cell * alpha_HEC_CD4T * alpha_CCL1921 * Treg_B...
    + gamma_LNB * c_LNB_cell * alpha_egress * Treg_LN - d_Treg * Treg_B - c_BPT_cell * (1 + alpha_TNFa_Treg * TNFa_PT / (K_TNFa_Treg + TNFa_PT)) * Treg_B;

% CTL
dydt(59) = -c_BLN_cell * alpha_HEC * alpha_CCL1921 * CTL_B...
    + gamma_LNB * c_LNB_cell * alpha_egress * CTL_LN - d_CTL * CTL_B - c_BPT_cell * (1 + alpha_TNFa_CTL * TNFa_PT / (K_TNFa_CTL + TNFa_PT)) * CTL_B;

% ASC
dydt(60) = -c_BLN_cell * alpha_HEC * alpha_CCL1921 * ASC_B...
    + gamma_LNB * c_LNB_cell * alpha_egress * ASC_LN - d_ASC * ASC_B - c_BPT_cell * (1 + alpha_TNFa_ASC * TNFa_PT / (K_TNFa_ASC + TNFa_PT)) * ASC_B;

% TNFa
dydt(61) = gamma_LNB * c_LNB_cyto * TNFa_LN + gamma_PTB * c_PTB_cyto * TNFa_PT - c_BLN_cyto * TNFa_B - c_BPT_cyto * TNFa_B - d_TNF * TNFa_B;

% IFNg
dydt(62) = gamma_LNB * c_LNB_cyto * IFNg_LN + gamma_PTB * c_PTB_cyto * IFNg_PT - c_BLN_cyto * IFNg_B - c_BPT_cyto * IFNg_B - d_IFNg * IFNg_B;

dydt(63) = 0;

% IL2
dydt(64) = gamma_LNB * c_LNB_cyto * IL2_LN - c_BLN_cyto * IL2_B + gamma_PTB * c_PTB_cyto * IL2_PT - c_BPT_cyto * IL2_B - d_IL2 * IL2_B;

% IL4
dydt(65) = gamma_LNB * c_LNB_cyto * IL4_LN - c_BLN_cyto * IL4_B + gamma_PTB * c_PTB_cyto * IL4_PT - c_BPT_cyto * IL4_B - d_IL4 * IL4_B;

% IL6
dydt(66) = gamma_LNB * c_LNB_cyto * IL6_LN - c_BLN_cyto * IL6_B + gamma_PTB * c_PTB_cyto * IL6_PT - c_BPT_cyto * IL6_B - d_IL6 * IL6_B;

% IL10
dydt(67) = (gamma_LNB * c_LNB_cyto * IL10_LN + gamma_PTB * c_PTB_cyto * IL10_PT - c_BLN_cyto * IL10_B - c_BPT_cyto * IL10_B - d_IL10 ...
    / (1 + (IL10_B / K_IL10d)) * IL10_B);

% Ab
dydt(68) = gamma_LNB * c_LNB_cyto * Ab_LN + gamma_PTB * c_PTB_cyto * Ab_PT - c_BLN_cyto * Ab_B - c_BPT_cyto * Ab_B - d_Ab * Ab_B;

% CA_B(TGF-β)
dydt(69) = gamma_PTB * c_PTB_cyto * CA_PT + gamma_LNB * c_LNB_cyto * CA_LN - c_BPT_cyto * CA_B - c_BLN_cyto * CA_B - d_CA * CA_B;

% D
dydt(70) = a_D * IL6_PT ^ 4 / (IL6_PT ^ 4 + K_D ^ 4) + a_P * Ag_PT / (Ag_PT + K_P ) - d_D * D;

%CA_PT(TGF-β)
dydt(71) = a_CA * mphi_PT - d_CA * CA_PT + s_CA + gamma_BPT * c_PTB_cyto * CA_B - c_PTB_cyto * CA_PT - c_PTLN_cyto * CA_PT;

%YIL10
dydt(72) = a_IL102 * D ^ 4 / (D ^ 4 + K_IL102 ^ 4) - d_IL102 * Y_IL10;

% equations of clock model
% mRNA_Per
dydt(73) = -1 * dm_per * Per + vmax_per * ( ...
    1 + fold_per * power(CLOCK_BMAL1 / Ka_per_cb, hill_per_cb)) / ( ...
    1 + power(CLOCK_BMAL1 / Ka_per_cb, hill_per_cb) * ...
    (1 + power(PER_CRY / Ki_per_pc, hill_per_pc)));

% mRNA_Cry
dydt(74) = -1 * dm_cry * Cry + 1 / ( ...
    1 + power(REV / Ki_cry_rev, hill_cry_rev)) * ( ...
    vmax_cry * ...
    (1 + fold_cry * power(CLOCK_BMAL1 / Ka_cry_cb, hill_cry_cb))) / ( ...
    1 + power(CLOCK_BMAL1 / Ka_cry_cb, hill_cry_cb) * ...
    (1 + power(PER_CRY / Ki_cry_pc, hill_cry_pc)));

% mRNA_Rev
dydt(75) = -1 * dm_rev * Rev + vmax_rev * ( ...
    1 + fold_rev * power(CLOCK_BMAL1 / Ka_rev_cb, hill_rev_cb)) / ( ...
    1 + power(CLOCK_BMAL1 / Ka_rev_cb, hill_rev_cb) * ...
    (1 + power(PER_CRY / Ki_rev_pc, hill_rev_pc)));

% mRNA_Ror
dydt(76) = -1 * dm_ror * Ror + vmax_ror * ( ...
    1 + fold_ror * power(CLOCK_BMAL1 / Ka_ror_cb, hill_ror_cb)) / ( ...
    1 + power(CLOCK_BMAL1 / Ka_ror_cb, hill_ror_cb) * ...
    (1 + power(PER_CRY / Ki_ror_pc, hill_ror_pc)));

% mRNA_Bmal
dydt(77) = -1 * dm_bmal * Bmal1 + vmax_bmal * K_P / (K_P + 0) * ( ...
    1 + fold_bmal * power(ROR / Ka_bmal_ror, hill_bmal_ror)) / ( ...
    1 + power(REV / Ki_bmal_rev, hill_bmal_rev) + ...
    power(ROR / Ka_bmal_ror, hill_bmal_ror));

% protein_PER
dydt(78) = control_PER * (-1 * dp_per * PER + kp_per * Per - ( ...
    kass_pc * PER * CRY - kdiss_pc * PER_CRY));

% protein_CRY
dydt(79) = control_CRY * (-1 * dp_cry * CRY + kp_cry * Cry - ...
    (kass_pc * PER * CRY - kdiss_pc * PER_CRY));

% dydt(77) = -1 * dp_rev * (1 + alpha_IL1b_REV * IL1b_PT .^ n_IL1b_REV / (IL1b_PT .^ n_IL1b_REV + K_IL1b_REV .^ n_IL1b_REV )) * REV + kp_rev * Rev;
% protein_REV-ERB
dydt(80) = control_REV * (-1 * dp_rev * REV + kp_rev * Rev);

% protein_ROR
dydt(81) = control_ROR * (-1 * dp_ror * ROR + kp_ror * Ror);

% protein_BMAL1
dydt(82) = control_BMAL1 * (-1 * dp_bmal * BMAL1 + kp_bmal * Bmal1 - ...
    (kass_cb * BMAL1 - kdiss_cb * CLOCK_BMAL1));

% protein_PER:CRY
dydt(83) = (kass_pc * PER * CRY - kdiss_pc * PER_CRY - d_pc * PER_CRY);

% protein_CLOCK:BMAL1
dydt(84) = (kass_cb * BMAL1 - kdiss_cb * CLOCK_BMAL1 - d_cb * CLOCK_BMAL1);

dydt = real(dydt);
if control_bistable_innate == 1
    dydt([3 : 5, 13 : 16, 20, 21, 23, 29, 42 : 44, 47, 48, 51, 61, 62, 66, 67, 69 : 72]) = 0;
end
end


function y = inflammation(DC, Mphi, IL6, TNF, IL10, TGF, Par)
%     y = (Par(1).*DC+Par(2).*Mphi+Par(3).*IL6+Par(4).*TNF).*(Par(5)./(IL10+1)+Par(6)./(TGF+1));
y = (Par(1) .* DC + Par(2) .* Mphi + Par(3) .* IL6 + Par(4) .* TNF);
end
function y = Ada_response(CTL, CD8m, Ab, Bm, IL2, IL4, Par)
%     y = (log10(Par(1)+1).*log10(CTL+1)+log10(Par(2)+1).*log10(CD8m+1)+log10(Par(3)+1).*log10(IL2+1)).*...
%         (log10(Par(4)+1).*log10(Ab+1)+log10(Par(5)+1).*log10(Bm+1)+log10(Par(6)+1).*log10(IL4+1));
%         y = (Par(1).*CTL+Par(2).*CD8m+Par(3).*IL2).*(Par(4).*Ab+Par(5).*Bm+Par(6).*IL4);
y = (Par(1) .* CTL + Par(2) .* CD8m + Par(3) .* IL2) + (Par(4) .* Ab + Par(5) .* Bm + Par(6) .* IL4);
end
