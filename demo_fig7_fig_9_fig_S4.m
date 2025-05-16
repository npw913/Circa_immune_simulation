function demo_fig7_fig_9_fig_S4
Par_decay = xlsread('decay_list.xlsx');
Par_produce = xlsread('produce_list.xlsx');
Par_transition = xlsread('transition_list.xlsx');
Par_halfmax = xlsread('half_maximum_list.xlsx');
Par_Others = xlsread('Others_list.xlsx');
Par_promotion = xlsread('promotion_list.xlsx');
Par_clock = load('par_clock.csv');
Par_Clock_immune = xlsread('circa_immune_list.xlsx');
Time_circadian_clock = xlsread('Timing.xlsx');
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

% Parameters of DC-CCL21 loop
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


% Phase diagram (DC vaccine)

% Time course of DC_PT&DC_LN&CCL21 (Generating the trajectories of ZT7 and ZT19)

T7_0 = 0 : dt : 7 / 24;
T19_0 = 0 : dt : (19/ 24);
T = 0 : dt : 4;

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
    Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,7), T, y7_0);
[~, Y19] = ode45(@(t, y)ODE_zero_circa_immune_response_fit_acute_version_T(t, y, Par_decay, Par_produce, ...
    Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,19), T, y19_0);
DCPT_ZT7 = Y7(:,3);
DCPT_ZT19 = Y19(:,3);
DCLN_ZT7 = Y7(:,23);
DCLN_ZT19 = Y19(:,23);
sCCL_ZT7 = Y7(:,50);
sCCL_ZT19 = Y19(:,50);
% fix point solution
        DC = 0:0.001:150;
        TEND=50:50:4000;
        fix_point_ensemble_DC=zeros(length(TEND),3);
        fix_point_ensemble_CCL=zeros(length(TEND),3);
        for ii = 1:length(TEND)
            tend=TEND(ii);
            by = Par_new(8) / Par_new(10); %8/10
            Kd = Par_new(9); %9
            bx =Par_new(1)*Par_new(2)*Par_new(5)*DCPT_ZT7(tend)/Par_new(6);  %1*2*Y19(4,tend)/6
            Kc = Par_new(3); 
            k=1;
            for j = 1:length(DC)-1
                 CCL_j =  by*DC(j).^2./(Kd^2+DC(j).^2);
                 CCL_jplus1 = by*DC(j+1).^2./(Kd^2+DC(j+1).^2);
                 a = DC(j) -  bx * CCL_j.^2./(Kc^2+CCL_j.^2);
                 b = DC(j+1) -  bx * CCL_jplus1.^2./(Kc^2+CCL_jplus1.^2);
                 if a*b<0
                     fix_point_ensemble_DC(ii,k)=(DC(j)+DC(j+1))/2;
                     fix_point_ensemble_CCL(ii,k)=(CCL_j+CCL_jplus1)/2;
                     k=k+1;
                 elseif a==0
                     fix_point_ensemble_DC(ii,k)=DC(j);
                     fix_point_ensemble_CCL(ii,k)=CCL_j;
                 elseif b==0
                     fix_point_ensemble_DC(ii,k)=DC(j+1);
                     fix_point_ensemble_CCL(ii,k)=CCL_jplus1;
                 else
                 end
            end

        end
f= figure("Name",'Fig7',NumberTitle="off");
TTT = tiledlayout(2, 2); 
TEND =[700 1400 1900 3000];
for jjj=1:4
    tend = TEND(jjj);
% tend = 1700;

set(gcf, 'Position', [100 100 800 600]); % 设置图形大小[宽度, 高度]
f.Color = [1 1 1];
%
nexttile(jjj)
% %
% % %     q=quiver(X_CCL,Y_DC,DyDt_xCCL,DyDt_yDC);
% % %     q.AutoScaleFactor=1;
        by = Par_new(8) / Par_new(10); %8/10
        Kd = Par_new(9); %9
        bx =Par_new(1)*Par_new(2)*Par_new(5)*DCPT_ZT7(tend)/Par_new(6);  %1*2*Y19(4,tend)/6
        Kc = Par_new(3); 
        DC_1 = 0:0.01:150;
%         ccl_1 = b0 + by*DC_1.^2./(Kd^2+DC_1.^2);
        ccl_1 = by*DC_1.^2./(Kd^2+DC_1.^2);
        ccl_2 = 0:0.01:800;
        DC_2 = bx * (ccl_2).^2./(Kc^2+(ccl_2).^2);
        p1=plot(ccl_1,DC_1,'LineStyle','--','Color','#800080','LineWidth',0.9);
        hold on
        p2=plot(ccl_2,DC_2,'LineStyle','--','Color','#008000','LineWidth',0.9);
        angle = -15;
        arrow_v = [1+0.02*cos(pi/180*angle) 1+0.02*sin(pi/180*angle)];
        arrow_nv = [1+0.02*cos(pi/180*(angle+180)) 1+0.02*sin(pi/180*(angle+180))];
        arrow_p90v = [1+0.02*cos(pi/180*(angle+45)) 1+0.02*sin(pi/180*(angle+45))];
        arrow_n90v = [1+0.02*cos(pi/180*(angle-135)) 1+0.02*sin(pi/180*(angle-135))];
    y0l = [fix_point_ensemble_CCL(tend/50,1)*arrow_nv(1),fix_point_ensemble_DC(tend/50,1)*arrow_nv(2),DCPT_ZT7(tend)];
    y0r = [fix_point_ensemble_CCL(tend/50,1)*arrow_v(1),fix_point_ensemble_DC(tend/50,1)*arrow_v(2),DCPT_ZT7(tend)];
    y0u = [fix_point_ensemble_CCL(tend/50,1)*arrow_p90v(1),fix_point_ensemble_DC(tend/50,1)*arrow_p90v(2),DCPT_ZT7(tend)];
    y0d = [fix_point_ensemble_CCL(tend/50,1)*arrow_n90v(1),fix_point_ensemble_DC(tend/50,1)*arrow_n90v(2),DCPT_ZT7(tend)];
    T  = 0:0.001:30;
    if tend<1900
    [~, Ytrl] = ode15s(@(t, y)ode_DCCCL_time_reversal(t, y, Par_new, Time_circadian_clock, 0), T, y0l);
    [~, Ytrr] = ode15s(@(t, y)ode_DCCCL_time_reversal(t, y, Par_new, Time_circadian_clock, 0), T, y0r);
    [~, Yntru] = ode15s(@(t, y)ode_DCCCL_no_time_reverse(t, y, Par_new, Time_circadian_clock, 0), T, y0u);
    [~, Yntrd] = ode15s(@(t, y)ode_DCCCL_no_time_reverse(t, y, Par_new, Time_circadian_clock, 0), T, y0d);
    CCL_manifold_l = Ytrl(:,1);
    DC_manifold_l = Ytrl(:,2);
    CCL_manifold_r = Ytrr(:,1);
    DC_manifold_r = Ytrr(:,2);
    CCL_manifold_u = Yntru(:,1);
    DC_manifold_u = Yntru(:,2);
    CCL_manifold_d = Yntrd(:,1);
    DC_manifold_d = Yntrd(:,2);
    i =0;
    while 1
        i =i+1;
%         if CCL_manifold_r(i)>800
%             DC_manifold_r(i) = DC_manifold_r(i)+(DC_manifold_r(i-1)-DC_manifold_r(i))*(800-CCL_manifold_r(i))/(CCL_manifold_r(i-1)-CCL_manifold_r(i));
%             CCL_manifold_r(i) = 800;
%             CCL_manifold_r((i+1):end)=[];
%             DC_manifold_r((i+1):end)=[];
%             break
%         end
        if DC_manifold_r(i)<0
            CCL_manifold_r(i) = CCL_manifold_r(i) - DC_manifold_r(i)*(CCL_manifold_r(i-1)-CCL_manifold_r(i))/(DC_manifold_r(i-1)-DC_manifold_r(i));
            DC_manifold_r(i) = 0;
            CCL_manifold_r((i+1):end)=[];
            DC_manifold_r((i+1):end)=[];
            break
        end
    end
    i=0;
    while 1
        i =i+1;
        if CCL_manifold_l(i)<0
            DC_manifold_l(i) = DC_manifold_l(i) - CCL_manifold_l(i)*(DC_manifold_l(i-1)-DC_manifold_l(i))/(CCL_manifold_l(i-1)-CCL_manifold_l(i));
            CCL_manifold_l(i) = 0;
            CCL_manifold_l((i+1):end)=[];
            DC_manifold_l((i+1):end)=[];
            break
        end
    end

    fill_CCL = [CCL_manifold_r(length(CCL_manifold_r):(-1):1);CCL_manifold_l;0];
    fill_DC = [DC_manifold_r(length(CCL_manifold_r):(-1):1);DC_manifold_l;0];
    i =1;
    while max(CCL_manifold_r)>800

        if CCL_manifold_r(i)>800
            DC_manifold_r(i) = DC_manifold_r(i)+(DC_manifold_r(i-1)-DC_manifold_r(i))*(800-CCL_manifold_r(i))/(CCL_manifold_r(i-1)-CCL_manifold_r(i));
            CCL_manifold_r(i) = 800;
            CCL_manifold_r((i+1):end)=[];
            DC_manifold_r((i+1):end)=[];
            break
        end
        i =i+1;
    end
    end
    if tend<1900
        if tend ==1400
    p3=arrowPlot(CCL_manifold_l(length(CCL_manifold_l):(-1):1),DC_manifold_l(length(CCL_manifold_l):(-1):1),...
        'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 700 0 12],'size',800);
    p4=arrowPlot(CCL_manifold_r(length(CCL_manifold_r):(-1):1),DC_manifold_r(length(CCL_manifold_r):(-1):1),...
    'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 700 0 12],'size',800);
    p5=arrowPlot(CCL_manifold_u,DC_manifold_u,...
        'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 700 0 12],'size',800);
    p6=arrowPlot(CCL_manifold_d,DC_manifold_d,...
    'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 700 0 12],'size',800);
    p7=fill(fill_CCL,fill_DC,[255 127 14]/255,'FaceAlpha',0.2,'EdgeColor','none','LineWidth',1.5);
        else
 
    
    p3=arrowPlot(CCL_manifold_l(length(CCL_manifold_l):(-1):1),DC_manifold_l(length(CCL_manifold_l):(-1):1),...
        'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 750 0 12],'size',800);
    p4=arrowPlot(CCL_manifold_r(length(CCL_manifold_r):(-1):1),DC_manifold_r(length(CCL_manifold_r):(-1):1),...
    'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 750 0 12],'size',800);
    p5=arrowPlot(CCL_manifold_u,DC_manifold_u,...
        'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 750 0 30],'size',800);
    p6=arrowPlot(CCL_manifold_d,DC_manifold_d,...
    'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 750 0 12],'size',800);
    p7=fill(fill_CCL,fill_DC,[255 127 14]/255,'FaceAlpha',0.2,'EdgeColor','none','LineWidth',1.5);
       end
    end
    p8=plot(fix_point_ensemble_CCL(tend/50,1),fix_point_ensemble_DC(tend/50,1),...
            'o', 'Color', [0.2 0.2 0.2],'MarkerFaceColor','w','MarkerSize',8,'LineWidth',2);
    p9=plot(fix_point_ensemble_CCL(tend/50,2),fix_point_ensemble_DC(tend/50,2),...
            'o', 'Color', [0.2 0.2 0.2], 'LineWidth', 3,'MarkerFaceColor',[0.2 0.2 0.2]);
    p10=plot(fix_point_ensemble_CCL(tend/50,3),fix_point_ensemble_DC(tend/50,3),...
            'o', 'Color', [0.2 0.2 0.2], 'LineWidth', 3,'MarkerFaceColor',[0.2 0.2 0.2]);
    p11=plot(sCCL_ZT7(1:tend), DCLN_ZT7(1:tend), '-', 'Color', '#1F77B4', 'LineWidth', 1.5);

    p12=plot(sCCL_ZT7(tend), DCLN_ZT7(tend),'o', 'Color', '#1F77B4', 'LineWidth', 3,'MarkerFaceColor','#1F77B4');
    p13=plot(sCCL_ZT19(1:tend), DCLN_ZT19(1:tend), '-', 'Color', '#D62728', 'LineWidth', 1.5);
    p14=plot(sCCL_ZT19(tend), DCLN_ZT19(tend),'o', 'Color', '#D62728', 'LineWidth', 3,'MarkerFaceColor','#D62728');
    if jjj~=1
    axis([0 800 0 12]);
%     axis square
    end
    if jjj==1
    xlim([0 800])
%     axis square
    symlog_deepseek('y',0.2,'Ylim',[0 60])  
    end
    hold off
    xlabel('Soluble CCL21(pg/ml)',FontSize=14);
    ylabel('DCs in dLNs(10^3/\mul)',FontSize=14);
    show_tend = tend/1000;
    title(['Time after Injection = ',num2str(show_tend,'%.2f'),' d'],FontSize=14);
%     if jjj ==3  
%     legend([p1,p2,p3,p7,p8,p9,p11,p13],{'Nullcline1(d[sCCL21]/dt=0)','Nullcline2(d[DC_{LN}]/dt=0)',...
%         'Saddle Manifold','Basin of Attraction','Saddle Point','Stable Steady State','ZT7','ZT19'},...
%         'location','southoutside',FontSize=12,Orientation='horizontal',NumColumns=4);
%     end
    set(gca, 'FontName', 'Arial','FontSize',12);
    ax = gca;
    ax.LineWidth = 1.5; % 加粗边框
    ax.TickLength = [0.012, 0.012]; % 前者设置y轴 tick的长度，后者设置x轴
    set(ax, 'TickDir', 'in')
%     box off
    hold off
%     F=getframe(gcf);
%     I=frame2im(F);
%     [I,map]=rgb2ind(I,256);
%     pic_num = tend/50;
%     if pic_num == 1
%         imwrite(I,map,'test2.gif','gif', 'Loopcount',inf,'DelayTime',0.2);
%     else
%         imwrite(I,map,'test2.gif','gif','WriteMode','append','DelayTime',0.2);
%     end
%     pic_num = pic_num + 1;
end

% Phase diagram (Nucleic-acid-adjuvanted vaccine)
% LPS/CpG converter
Par_Others(35)=1;
% Time course of DC_PT&DC_LN&CCL21 (Generating the trajectories of ZT7 and ZT19)

T7_0 = 0 : dt : 7 / 24;
T19_0 = 0 : dt : (19/ 24);
T = 0 : dt : 4;

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
    Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,7), T, y7_0);
[~, Y19] = ode45(@(t, y)ODE_zero_circa_immune_response_fit_acute_version_T(t, y, Par_decay, Par_produce, ...
    Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,19), T, y19_0);
DCPT_ZT7 = Y7(:,3);
DCPT_ZT19 = Y19(:,3);
DCLN_ZT7 = Y7(:,23);
DCLN_ZT19 = Y19(:,23);
sCCL_ZT7 = Y7(:,50);
sCCL_ZT19 = Y19(:,50);
% fix point solution
        DC = 0:0.001:150;
        TEND=50:50:4000;
        fix_point_ensemble_DC=zeros(length(TEND),3);
        fix_point_ensemble_CCL=zeros(length(TEND),3);
        for ii = 1:length(TEND)
            tend=TEND(ii);
            by = Par_new(8) / Par_new(10); %8/10
            Kd = Par_new(9); %9
            bx =Par_new(1)*Par_new(2)*Par_new(5)*DCPT_ZT19(tend)/Par_new(6);  %1*2*Y19(4,tend)/6
            Kc = Par_new(3); 
            k=1;
            for j = 1:length(DC)-1
                 CCL_j =  by*DC(j).^2./(Kd^2+DC(j).^2);
                 CCL_jplus1 = by*DC(j+1).^2./(Kd^2+DC(j+1).^2);
                 a = DC(j) -  bx * CCL_j.^2./(Kc^2+CCL_j.^2);
                 b = DC(j+1) -  bx * CCL_jplus1.^2./(Kc^2+CCL_jplus1.^2);
                 if a*b<0
                     fix_point_ensemble_DC(ii,k)=(DC(j)+DC(j+1))/2;
                     fix_point_ensemble_CCL(ii,k)=(CCL_j+CCL_jplus1)/2;
                     k=k+1;
                 elseif a==0
                     fix_point_ensemble_DC(ii,k)=DC(j);
                     fix_point_ensemble_CCL(ii,k)=CCL_j;
                 elseif b==0
                     fix_point_ensemble_DC(ii,k)=DC(j+1);
                     fix_point_ensemble_CCL(ii,k)=CCL_jplus1;
                 else
                 end
            end

        end

f= figure("Name",'Fig9',NumberTitle="off");
f.Color = [1 1 1];
TTT = tiledlayout(2, 2); 
TEND =[700 1500 1900 3000];
for jjj=1:4
    tend = TEND(jjj);
% tend = 1700;

set(gcf, 'Position', [100 100 800 600]); % 设置图形大小[宽度, 高度]

%
nexttile(jjj)
% %
% % %     q=quiver(X_CCL,Y_DC,DyDt_xCCL,DyDt_yDC);
% % %     q.AutoScaleFactor=1;
        by = Par_new(8) / Par_new(10); %8/10
        Kd = Par_new(9); %9
        bx =Par_new(1)*Par_new(2)*Par_new(5)*DCPT_ZT19(tend)/Par_new(6);  %1*2*Y19(4,tend)/6
        Kc = Par_new(3); 
        DC_1 = 0:0.01:150;
%         ccl_1 = b0 + by*DC_1.^2./(Kd^2+DC_1.^2);
        ccl_1 = by*DC_1.^2./(Kd^2+DC_1.^2);
        ccl_2 = 0:0.01:800;
        DC_2 = bx * (ccl_2).^2./(Kc^2+(ccl_2).^2);
        p1=plot(ccl_1,DC_1,'LineStyle','--','Color','#800080','LineWidth',0.9);
        hold on
        p2=plot(ccl_2,DC_2,'LineStyle','--','Color','#008000','LineWidth',0.9);
        angle = -15;
        arrow_v = [1+0.02*cos(pi/180*angle) 1+0.02*sin(pi/180*angle)];
        arrow_nv = [1+0.02*cos(pi/180*(angle+180)) 1+0.02*sin(pi/180*(angle+180))];
        arrow_p90v = [1+0.02*cos(pi/180*(angle+45)) 1+0.02*sin(pi/180*(angle+45))];
        arrow_n90v = [1+0.02*cos(pi/180*(angle-135)) 1+0.02*sin(pi/180*(angle-135))];
    y0l = [fix_point_ensemble_CCL(tend/50,1)*arrow_nv(1),fix_point_ensemble_DC(tend/50,1)*arrow_nv(2),DCPT_ZT19(tend)];
    y0r = [fix_point_ensemble_CCL(tend/50,1)*arrow_v(1),fix_point_ensemble_DC(tend/50,1)*arrow_v(2),DCPT_ZT19(tend)];
    y0u = [fix_point_ensemble_CCL(tend/50,1)*arrow_p90v(1),fix_point_ensemble_DC(tend/50,1)*arrow_p90v(2),DCPT_ZT19(tend)];
    y0d = [fix_point_ensemble_CCL(tend/50,1)*arrow_n90v(1),fix_point_ensemble_DC(tend/50,1)*arrow_n90v(2),DCPT_ZT19(tend)];
    T  = 0:0.001:30;
    if tend<1900
    [~, Ytrl] = ode15s(@(t, y)ode_DCCCL_time_reversal(t, y, Par_new, Time_circadian_clock, 0), T, y0l);
    [~, Ytrr] = ode15s(@(t, y)ode_DCCCL_time_reversal(t, y, Par_new, Time_circadian_clock, 0), T, y0r);
    [~, Yntru] = ode15s(@(t, y)ode_DCCCL_no_time_reverse(t, y, Par_new, Time_circadian_clock, 0), T, y0u);
    [~, Yntrd] = ode15s(@(t, y)ode_DCCCL_no_time_reverse(t, y, Par_new, Time_circadian_clock, 0), T, y0d);
    CCL_manifold_l = Ytrl(:,1);
    DC_manifold_l = Ytrl(:,2);
    CCL_manifold_r = Ytrr(:,1);
    DC_manifold_r = Ytrr(:,2);
    CCL_manifold_u = Yntru(:,1);
    DC_manifold_u = Yntru(:,2);
    CCL_manifold_d = Yntrd(:,1);
    DC_manifold_d = Yntrd(:,2);
    i =0;
    while 1
        i =i+1;
%         if CCL_manifold_r(i)>800
%             DC_manifold_r(i) = DC_manifold_r(i)+(DC_manifold_r(i-1)-DC_manifold_r(i))*(800-CCL_manifold_r(i))/(CCL_manifold_r(i-1)-CCL_manifold_r(i));
%             CCL_manifold_r(i) = 800;
%             CCL_manifold_r((i+1):end)=[];
%             DC_manifold_r((i+1):end)=[];
%             break
%         end
        if DC_manifold_r(i)<0
            CCL_manifold_r(i) = CCL_manifold_r(i) - DC_manifold_r(i)*(CCL_manifold_r(i-1)-CCL_manifold_r(i))/(DC_manifold_r(i-1)-DC_manifold_r(i));
            DC_manifold_r(i) = 0;
            CCL_manifold_r((i+1):end)=[];
            DC_manifold_r((i+1):end)=[];
            break
        end
    end
    i=0;
    while 1
        i =i+1;
        if CCL_manifold_l(i)<0
            DC_manifold_l(i) = DC_manifold_l(i) - CCL_manifold_l(i)*(DC_manifold_l(i-1)-DC_manifold_l(i))/(CCL_manifold_l(i-1)-CCL_manifold_l(i));
            CCL_manifold_l(i) = 0;
            CCL_manifold_l((i+1):end)=[];
            DC_manifold_l((i+1):end)=[];
            break
        end
    end

    fill_CCL = [CCL_manifold_r(length(CCL_manifold_r):(-1):1);CCL_manifold_l;0];
    fill_DC = [DC_manifold_r(length(CCL_manifold_r):(-1):1);DC_manifold_l;0];
    i =1;
    while max(CCL_manifold_r)>800

        if CCL_manifold_r(i)>800
            DC_manifold_r(i) = DC_manifold_r(i)+(DC_manifold_r(i-1)-DC_manifold_r(i))*(800-CCL_manifold_r(i))/(CCL_manifold_r(i-1)-CCL_manifold_r(i));
            CCL_manifold_r(i) = 800;
            CCL_manifold_r((i+1):end)=[];
            DC_manifold_r((i+1):end)=[];
            break
        end
        i =i+1;
    end
    end
    if tend<1900
        if tend ==1500
    p3=arrowPlot(CCL_manifold_l(length(CCL_manifold_l):(-1):1),DC_manifold_l(length(CCL_manifold_l):(-1):1),...
        'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 700 0 12],'size',800);
    p4=arrowPlot(CCL_manifold_r(length(CCL_manifold_r):(-1):1),DC_manifold_r(length(CCL_manifold_r):(-1):1),...
    'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 700 0 12],'size',800);
    p5=arrowPlot(CCL_manifold_u,DC_manifold_u,...
        'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 700 0 12],'size',800);
    p6=arrowPlot(CCL_manifold_d,DC_manifold_d,...
    'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 700 0 12],'size',800);
    p7=fill(fill_CCL,fill_DC,[255 127 14]/255,'FaceAlpha',0.2,'EdgeColor','none','LineWidth',1.5);
        else
 
    
    p3=arrowPlot(CCL_manifold_l(length(CCL_manifold_l):(-1):1),DC_manifold_l(length(CCL_manifold_l):(-1):1),...
        'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 750 0 12],'size',800);
    p4=arrowPlot(CCL_manifold_r(length(CCL_manifold_r):(-1):1),DC_manifold_r(length(CCL_manifold_r):(-1):1),...
    'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 750 0 12],'size',800);
    p5=arrowPlot(CCL_manifold_u,DC_manifold_u,...
        'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 750 0 30],'size',800);
    p6=arrowPlot(CCL_manifold_d,DC_manifold_d,...
    'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 750 0 12],'size',800);
    p7=fill(fill_CCL,fill_DC,[255 127 14]/255,'FaceAlpha',0.2,'EdgeColor','none','LineWidth',1.5);
       end
    end
    p8=plot(fix_point_ensemble_CCL(tend/50,1),fix_point_ensemble_DC(tend/50,1),...
            'o', 'Color', [0.2 0.2 0.2],'MarkerFaceColor','w','MarkerSize',8,'LineWidth',2);
    p9=plot(fix_point_ensemble_CCL(tend/50,2),fix_point_ensemble_DC(tend/50,2),...
            'o', 'Color', [0.2 0.2 0.2], 'LineWidth', 3,'MarkerFaceColor',[0.2 0.2 0.2]);
    p10=plot(fix_point_ensemble_CCL(tend/50,3),fix_point_ensemble_DC(tend/50,3),...
            'o', 'Color', [0.2 0.2 0.2], 'LineWidth', 3,'MarkerFaceColor',[0.2 0.2 0.2]);
    p11=plot(sCCL_ZT7(1:tend), DCLN_ZT7(1:tend), '-', 'Color', '#1F77B4', 'LineWidth', 1.5);

    p12=plot(sCCL_ZT7(tend), DCLN_ZT7(tend),'o', 'Color', '#1F77B4', 'LineWidth', 3,'MarkerFaceColor','#1F77B4');
    p13=plot(sCCL_ZT19(1:tend), DCLN_ZT19(1:tend), '-', 'Color', '#D62728', 'LineWidth', 1.5);
    p14=plot(sCCL_ZT19(tend), DCLN_ZT19(tend),'o', 'Color', '#D62728', 'LineWidth', 3,'MarkerFaceColor','#D62728');
    if jjj~=1
    axis([0 800 0 12]);
%     axis square
    end
    if jjj==1
    xlim([0 800])
%     axis square
    symlog_deepseek('y',0.2,'Ylim',[0 60])  
    end
    hold off
    xlabel('Soluble CCL21(pg/ml)',FontSize=14);
    ylabel('DCs in dLNs(10^3/\mul)',FontSize=14);
    show_tend = tend/1000;
    title(['Time after Injection = ',num2str(show_tend,'%.2f'),' d'],FontSize=14);
%     if jjj ==3  
%     legend([p1,p2,p3,p7,p8,p9,p11,p13],{'Nullcline1(d[sCCL21]/dt=0)','Nullcline2(d[DC_{LN}]/dt=0)',...
%         'Saddle Manifold','Basin of Attraction','Saddle Point','Stable Steady State','ZT7','ZT19'},...
%         'location','southoutside',FontSize=12,Orientation='horizontal',NumColumns=4);
%     end
    set(gca, 'FontName', 'Arial','FontSize',12);
    ax = gca;
    ax.LineWidth = 1.5; % 加粗边框
    ax.TickLength = [0.012, 0.012]; % 前者设置y轴 tick的长度，后者设置x轴
    set(ax, 'TickDir', 'in')
%     box off
    hold off
%     F=getframe(gcf);
%     I=frame2im(F);
%     [I,map]=rgb2ind(I,256);
%     pic_num = tend/50;
%     if pic_num == 1
%         imwrite(I,map,'test2.gif','gif', 'Loopcount',inf,'DelayTime',0.2);
%     else
%         imwrite(I,map,'test2.gif','gif','WriteMode','append','DelayTime',0.2);
%     end
%     pic_num = pic_num + 1;
end

% Phase diagram (Endotoxin-adjuvanted vaccine)
Par_Others(35)=0;
% Time course of DC_PT&DC_LN&CCL21 (Generating the trajectories of ZT7 and ZT19)

T7_0 = 0 : dt : 7 / 24;
T19_0 = 0 : dt : (19/ 24);
T = 0 : dt : 4;

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
    Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,7), T, y7_0);
[~, Y19] = ode45(@(t, y)ODE_zero_circa_immune_response_fit_acute_version_T(t, y, Par_decay, Par_produce, ...
    Par_transition, Par_halfmax, Par_Others, Par_promotion, Par_clock, Par_Clock_immune,Time_circadian_clock,19), T, y19_0);
DCPT_ZT7 = Y7(:,3);
DCPT_ZT19 = Y19(:,3);
DCLN_ZT7 = Y7(:,23);
DCLN_ZT19 = Y19(:,23);
sCCL_ZT7 = Y7(:,50);
sCCL_ZT19 = Y19(:,50);
% fix point solution
        DC = 0:0.001:150;
        TEND=50:50:4000;
        fix_point_ensemble_DC=zeros(length(TEND),3);
        fix_point_ensemble_CCL=zeros(length(TEND),3);
        for ii = 1:length(TEND)
            tend=TEND(ii);
            by = Par_new(8) / Par_new(10); %8/10
            Kd = Par_new(9); %9
            bx =Par_new(1)*Par_new(2)*Par_new(5)*DCPT_ZT7(tend)/Par_new(6);  %1*2*Y19(4,tend)/6
            Kc = Par_new(3); 
            k=1;
            for j = 1:length(DC)-1
                 CCL_j =  by*DC(j).^2./(Kd^2+DC(j).^2);
                 CCL_jplus1 = by*DC(j+1).^2./(Kd^2+DC(j+1).^2);
                 a = DC(j) -  bx * CCL_j.^2./(Kc^2+CCL_j.^2);
                 b = DC(j+1) -  bx * CCL_jplus1.^2./(Kc^2+CCL_jplus1.^2);
                 if a*b<0
                     fix_point_ensemble_DC(ii,k)=(DC(j)+DC(j+1))/2;
                     fix_point_ensemble_CCL(ii,k)=(CCL_j+CCL_jplus1)/2;
                     k=k+1;
                 elseif a==0
                     fix_point_ensemble_DC(ii,k)=DC(j);
                     fix_point_ensemble_CCL(ii,k)=CCL_j;
                 elseif b==0
                     fix_point_ensemble_DC(ii,k)=DC(j+1);
                     fix_point_ensemble_CCL(ii,k)=CCL_jplus1;
                 else
                 end
            end

        end
f= figure("Name",'FigS4',NumberTitle="off");
TTT = tiledlayout(2, 2); 
TEND =[700 1400 1900 3000];
for jjj=1:4
    tend = TEND(jjj);
% tend = 1700;

set(gcf, 'Position', [100 100 800 600]); % 设置图形大小[宽度, 高度]
f.Color = [1 1 1];
%
nexttile(jjj)
% %
% % %     q=quiver(X_CCL,Y_DC,DyDt_xCCL,DyDt_yDC);
% % %     q.AutoScaleFactor=1;
        by = Par_new(8) / Par_new(10); %8/10
        Kd = Par_new(9); %9
        bx =Par_new(1)*Par_new(2)*Par_new(5)*DCPT_ZT7(tend)/Par_new(6);  %1*2*Y19(4,tend)/6
        Kc = Par_new(3); 
        DC_1 = 0:0.01:150;
%         ccl_1 = b0 + by*DC_1.^2./(Kd^2+DC_1.^2);
        ccl_1 = by*DC_1.^2./(Kd^2+DC_1.^2);
        ccl_2 = 0:0.01:800;
        DC_2 = bx * (ccl_2).^2./(Kc^2+(ccl_2).^2);
        p1=plot(ccl_1,DC_1,'LineStyle','--','Color','#800080','LineWidth',0.9);
        hold on
        p2=plot(ccl_2,DC_2,'LineStyle','--','Color','#008000','LineWidth',0.9);
        angle = -15;
        arrow_v = [1+0.02*cos(pi/180*angle) 1+0.02*sin(pi/180*angle)];
        arrow_nv = [1+0.02*cos(pi/180*(angle+180)) 1+0.02*sin(pi/180*(angle+180))];
        arrow_p90v = [1+0.02*cos(pi/180*(angle+45)) 1+0.02*sin(pi/180*(angle+45))];
        arrow_n90v = [1+0.02*cos(pi/180*(angle-135)) 1+0.02*sin(pi/180*(angle-135))];
    y0l = [fix_point_ensemble_CCL(tend/50,1)*arrow_nv(1),fix_point_ensemble_DC(tend/50,1)*arrow_nv(2),DCPT_ZT7(tend)];
    y0r = [fix_point_ensemble_CCL(tend/50,1)*arrow_v(1),fix_point_ensemble_DC(tend/50,1)*arrow_v(2),DCPT_ZT7(tend)];
    y0u = [fix_point_ensemble_CCL(tend/50,1)*arrow_p90v(1),fix_point_ensemble_DC(tend/50,1)*arrow_p90v(2),DCPT_ZT7(tend)];
    y0d = [fix_point_ensemble_CCL(tend/50,1)*arrow_n90v(1),fix_point_ensemble_DC(tend/50,1)*arrow_n90v(2),DCPT_ZT7(tend)];
    T  = 0:0.001:30;
    if tend<1900
    [~, Ytrl] = ode15s(@(t, y)ode_DCCCL_time_reversal(t, y, Par_new, Time_circadian_clock, 0), T, y0l);
    [~, Ytrr] = ode15s(@(t, y)ode_DCCCL_time_reversal(t, y, Par_new, Time_circadian_clock, 0), T, y0r);
    [~, Yntru] = ode15s(@(t, y)ode_DCCCL_no_time_reverse(t, y, Par_new, Time_circadian_clock, 0), T, y0u);
    [~, Yntrd] = ode15s(@(t, y)ode_DCCCL_no_time_reverse(t, y, Par_new, Time_circadian_clock, 0), T, y0d);
    CCL_manifold_l = Ytrl(:,1);
    DC_manifold_l = Ytrl(:,2);
    CCL_manifold_r = Ytrr(:,1);
    DC_manifold_r = Ytrr(:,2);
    CCL_manifold_u = Yntru(:,1);
    DC_manifold_u = Yntru(:,2);
    CCL_manifold_d = Yntrd(:,1);
    DC_manifold_d = Yntrd(:,2);
    i =0;
    while 1
        i =i+1;
%         if CCL_manifold_r(i)>800
%             DC_manifold_r(i) = DC_manifold_r(i)+(DC_manifold_r(i-1)-DC_manifold_r(i))*(800-CCL_manifold_r(i))/(CCL_manifold_r(i-1)-CCL_manifold_r(i));
%             CCL_manifold_r(i) = 800;
%             CCL_manifold_r((i+1):end)=[];
%             DC_manifold_r((i+1):end)=[];
%             break
%         end
        if DC_manifold_r(i)<0
            CCL_manifold_r(i) = CCL_manifold_r(i) - DC_manifold_r(i)*(CCL_manifold_r(i-1)-CCL_manifold_r(i))/(DC_manifold_r(i-1)-DC_manifold_r(i));
            DC_manifold_r(i) = 0;
            CCL_manifold_r((i+1):end)=[];
            DC_manifold_r((i+1):end)=[];
            break
        end
    end
    i=0;
    while 1
        i =i+1;
        if CCL_manifold_l(i)<0
            DC_manifold_l(i) = DC_manifold_l(i) - CCL_manifold_l(i)*(DC_manifold_l(i-1)-DC_manifold_l(i))/(CCL_manifold_l(i-1)-CCL_manifold_l(i));
            CCL_manifold_l(i) = 0;
            CCL_manifold_l((i+1):end)=[];
            DC_manifold_l((i+1):end)=[];
            break
        end
    end

    fill_CCL = [CCL_manifold_r(length(CCL_manifold_r):(-1):1);CCL_manifold_l;0];
    fill_DC = [DC_manifold_r(length(CCL_manifold_r):(-1):1);DC_manifold_l;0];
    i =1;
    while max(CCL_manifold_r)>800

        if CCL_manifold_r(i)>800
            DC_manifold_r(i) = DC_manifold_r(i)+(DC_manifold_r(i-1)-DC_manifold_r(i))*(800-CCL_manifold_r(i))/(CCL_manifold_r(i-1)-CCL_manifold_r(i));
            CCL_manifold_r(i) = 800;
            CCL_manifold_r((i+1):end)=[];
            DC_manifold_r((i+1):end)=[];
            break
        end
        i =i+1;
    end
    end
    if tend<1900
        if tend ==1400
    p3=arrowPlot(CCL_manifold_l(length(CCL_manifold_l):(-1):1),DC_manifold_l(length(CCL_manifold_l):(-1):1),...
        'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 700 0 12],'size',800);
    p4=arrowPlot(CCL_manifold_r(length(CCL_manifold_r):(-1):1),DC_manifold_r(length(CCL_manifold_r):(-1):1),...
    'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 700 0 12],'size',800);
    p5=arrowPlot(CCL_manifold_u,DC_manifold_u,...
        'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 700 0 12],'size',800);
    p6=arrowPlot(CCL_manifold_d,DC_manifold_d,...
    'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 700 0 12],'size',800);
    p7=fill(fill_CCL,fill_DC,[255 127 14]/255,'FaceAlpha',0.2,'EdgeColor','none','LineWidth',1.5);
        else
 
    
    p3=arrowPlot(CCL_manifold_l(length(CCL_manifold_l):(-1):1),DC_manifold_l(length(CCL_manifold_l):(-1):1),...
        'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 750 0 12],'size',800);
    p4=arrowPlot(CCL_manifold_r(length(CCL_manifold_r):(-1):1),DC_manifold_r(length(CCL_manifold_r):(-1):1),...
    'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 750 0 12],'size',800);
    p5=arrowPlot(CCL_manifold_u,DC_manifold_u,...
        'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 750 0 30],'size',800);
    p6=arrowPlot(CCL_manifold_d,DC_manifold_d,...
    'number', 1, 'color', [0.5 0.5 0.5], 'LineWidth', 0.9, 'scale', 1,'limit',[0 750 0 12],'size',800);
    p7=fill(fill_CCL,fill_DC,[255 127 14]/255,'FaceAlpha',0.2,'EdgeColor','none','LineWidth',1.5);
       end
    end
    p8=plot(fix_point_ensemble_CCL(tend/50,1),fix_point_ensemble_DC(tend/50,1),...
            'o', 'Color', [0.2 0.2 0.2],'MarkerFaceColor','w','MarkerSize',8,'LineWidth',2);
    p9=plot(fix_point_ensemble_CCL(tend/50,2),fix_point_ensemble_DC(tend/50,2),...
            'o', 'Color', [0.2 0.2 0.2], 'LineWidth', 3,'MarkerFaceColor',[0.2 0.2 0.2]);
    p10=plot(fix_point_ensemble_CCL(tend/50,3),fix_point_ensemble_DC(tend/50,3),...
            'o', 'Color', [0.2 0.2 0.2], 'LineWidth', 3,'MarkerFaceColor',[0.2 0.2 0.2]);
    p11=plot(sCCL_ZT7(1:tend), DCLN_ZT7(1:tend), '-', 'Color', '#1F77B4', 'LineWidth', 1.5);

    p12=plot(sCCL_ZT7(tend), DCLN_ZT7(tend),'o', 'Color', '#1F77B4', 'LineWidth', 3,'MarkerFaceColor','#1F77B4');
    p13=plot(sCCL_ZT19(1:tend), DCLN_ZT19(1:tend), '-', 'Color', '#D62728', 'LineWidth', 1.5);
    p14=plot(sCCL_ZT19(tend), DCLN_ZT19(tend),'o', 'Color', '#D62728', 'LineWidth', 3,'MarkerFaceColor','#D62728');
    if jjj~=1
    axis([0 800 0 12]);
%     axis square
    end
    if jjj==1
    xlim([0 800])
%     axis square
    symlog_deepseek('y',0.2,'Ylim',[0 60])  
    end
    hold off
    xlabel('Soluble CCL21(pg/ml)',FontSize=14);
    ylabel('DCs in dLNs(10^3/\mul)',FontSize=14);
    show_tend = tend/1000;
    title(['Time after Injection = ',num2str(show_tend,'%.2f'),' d'],FontSize=14);
%     if jjj ==3  
%     legend([p1,p2,p3,p7,p8,p9,p11,p13],{'Nullcline1(d[sCCL21]/dt=0)','Nullcline2(d[DC_{LN}]/dt=0)',...
%         'Saddle Manifold','Basin of Attraction','Saddle Point','Stable Steady State','ZT7','ZT19'},...
%         'location','southoutside',FontSize=12,Orientation='horizontal',NumColumns=4);
%     end
    set(gca, 'FontName', 'Arial','FontSize',12);
    ax = gca;
    ax.LineWidth = 1.5; % 加粗边框
    ax.TickLength = [0.012, 0.012]; % 前者设置y轴 tick的长度，后者设置x轴
    set(ax, 'TickDir', 'in')
%     box off
    hold off
%     F=getframe(gcf);
%     I=frame2im(F);
%     [I,map]=rgb2ind(I,256);
%     pic_num = tend/50;
%     if pic_num == 1
%         imwrite(I,map,'test2.gif','gif', 'Loopcount',inf,'DelayTime',0.2);
%     else
%         imwrite(I,map,'test2.gif','gif','WriteMode','append','DelayTime',0.2);
%     end
%     pic_num = pic_num + 1;
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
    CCL21_LN, S1PR1, TNFa_LN, IFNg_LN, CA_LN, IL2_LN, IL4_LN, IL6_LN, IL10_LN, Ab_LN, ...       % 40~49
    CCL19_LN, Neu_B, CD4Tn_B, CD8Tn_B, Bn_B, Th_B, Th2_B, Th17_B, Treg_B, CTL_B, ASC_B, ...          % 50~60
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

%引入circa_周期函数的目的：为了实现特异性敲除
%稳态时，ZT_initial统一取0
%响应态时，ZT_initial取注射时刻的ZT
% 要保证t大于等于0，ZT0最好取24，避免报错
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
% 抗原
dydt(1) = - d_Ag * Ag_PT;
% 免疫佐剂
dydt(2) = - d_Ad * Ad_PT;
% * 1 / (1 + (ROR_t / K_TNFROR))
%DC

if control_CpG == 1
    R_DC = alpha_DCP * Ag_PT * 1 / (1 + (CA_PT / K_DCCA)) * 1 / (1 + (IL10_PT / K_DCIL10)) * 2 / (1 + (CLOCK_BMAL1_t / K_BC_APC) ^ 3) * (1 + alpha_DCTNF * TNFa_PT / (TNFa_PT + K_DCTNF )...
        + alpha_DCIFN * IFNg_PT / (IFNg_PT + K_DCIFN ));
    dydt(3) = a_DC * R_DC / (K_DC + R_DC) - c_PTLN_cell * b * (CCL1921_DC_LN ^ n_CCL1921_APCm) / (K_CCL1921_APCm ^ n_CCL1921_APCm + (CCL1921_DC_LN) ^ n_CCL1921_APCm) * DC_PT  ...
        - d_DC * DC_PT;
end
if control_CpG == 0
    %     if Ag_PT>1e-50
    %        R_DC = alpha_DCP *12*exp(-16*t) * 1 / (1 + (CA_PT / K_DCCA)) * 1 / (1 + (IL10_PT / K_DCIL10)) * (1 + alpha_DCTNF * TNFa_PT / (TNFa_PT + K_DCTNF )...
    %         + alpha_DCIFN * IFNg_PT / (IFNg_PT + K_DCIFN ));
    %        dydt(3) = a_DC * R_DC / (K_DC + R_DC) - c_PTLN_cell * b * (CCL1921_DC_LN ^ n_CCL1921_APCm) / (K_CCL1921_APCm ^ n_CCL1921_APCm + (CCL1921_DC_LN) ^ n_CCL1921_APCm) * DC_PT  ...
    %         - d_DC * DC_PT;
    %     else
    %         dydt(3) =- c_PTLN_cell * b * (CCL1921_DC_LN ^ n_CCL1921_APCm) / (K_CCL1921_APCm ^ n_CCL1921_APCm + (CCL1921_DC_LN) ^ n_CCL1921_APCm) * DC_PT- d_DC * DC_PT;
    %     end
    R_DC = alpha_DCP * Ag_PT * 1 / (1 + (CA_PT / K_DCCA)) * 1 / (1 + (IL10_PT / K_DCIL10)) * (1 + alpha_DCTNF * TNFa_PT / (TNFa_PT + K_DCTNF )...
        + alpha_DCIFN * IFNg_PT / (IFNg_PT + K_DCIFN ));
    % % % %     R_DC = alpha_DCP * Ag_PT;
    dydt(3) = a_DC * R_DC / (K_DC + R_DC) - c_PTLN_cell * b * (CCL1921_DC_LN ^ n_CCL1921_APCm) / (K_CCL1921_APCm ^ n_CCL1921_APCm + (CCL1921_DC_LN) ^ n_CCL1921_APCm) * DC_PT  ...
        - d_DC * DC_PT;
    % %     dydt(3) = 0;
end
% Ag, no Ad.
% dydt(3) = gamma_BPT * c_BPT_APCi * (1 + alpha_TNFa_APCi * TNFa_PT / (K_TNFa_APCi + TNFa_PT)) * CCL2_PT ^ 2 / (K_CCL2_APCi ^ 2 + CCL2_PT ^ 2) * APCi_B - ...
%     d_APCi * APCi_PT - a_Ag_APCm * (1 + alpha_IFNg_APCi * IFNg_PT / (K_IFNg_APCi + IFNg_PT) ) * K_IL10_APCi / (K_IL10_APCi + IL10_PT) * ...
%     K_REV_APCi ^ 2 / (K_REV_APCi ^ 2 + REV ^ 2) * APCi_PT * Ag_PT;

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
% dydt(5) = 0;
%Treg
dydt(6) = gamma_BPT * c_BPT_cell * (1 + alpha_TNFa_Treg * TNFa_PT / (K_TNFa_Treg + TNFa_PT)) * Treg_B - c_PTLN_cell * Treg_PT - d_Treg * Treg_PT;

%Th
dydt(7) = gamma_BPT * c_BPT_cell * (1 + alpha_TNFa_Th * TNFa_PT / (K_TNFa_Th + TNFa_PT)) * Th_B - c_PTLN_cell * Th_PT - d_Th * Th_PT - d_Treg_T * Th_PT * Treg_PT;

% dydt(7) = gamma_BPT * c_BPT_Th1 * (1 + alpha_TNFa_Th1 * TNFa_PT / (K_TNFa_Th1 + TNFa_PT)) * Th1_B - c_PTLN_Th1 * Th1_PT - d_Th * Th1_PT - d_Treg_T * Th1_PT * Treg_PT;

% mRNA_TNFa_DC_PT
dydt(8) = 0;


dydt(9) = 0;

% CTL
dydt(10) = gamma_BPT * c_BPT_cell * (1 + alpha_TNFa_CTL * TNFa_PT / (K_TNFa_CTL + TNFa_PT)) * CTL_B - c_PTLN_cell * CTL_PT - d_CTL * CTL_PT - d_Treg_T * CTL_PT * Treg_PT;

% ASC
dydt(11) = gamma_BPT * c_BPT_cell * (1 + alpha_TNFa_ASC * TNFa_PT / (K_TNFa_ASC + TNFa_PT)) * ASC_B - c_PTLN_cell * ASC_PT - d_ASC * ASC_PT;

% CCL2不需要
dydt(12) = 0;

% CXCL5
dydt(13) = a_0_CXCL5 + a_Neu_CXCL5 * K_IL10_CXCL5 / ( IL10_PT + K_IL10_CXCL5) * Neu_PT - d_CXCL5 * CXCL5_PT;
%  dydt(13) = 0;
% CXCL8
dydt(14) = a_APCm_CXCL8 * K_IL10_CXCL8 / ( IL10_PT + K_IL10_CXCL8) * mphi_PT - d_CXCL8 * CXCL8_PT;
% dydt(14) = 0;
% dydt(15) = a_APCm_TNFa * (1 + alpha_IL1b_TNFa * IL1b_PT / ( IL1b_PT + K_IL1b_TNFa)) * K_IL10_TNFa / ( IL10_PT + K_IL10_TNFa)...
%     * K_ROR_TNFa / (K_ROR_TNFa + ROR) * APCm_PT + ...
%     a_Th1_TNFa * Th1_PT + a_CTL_TNFa * CTL_PT + gamma_BPT * c_BPT_TNFa * TNFa_B - c_PTLN_TNFa * TNFa_PT - d_TNFa * TNFa_PT;
% dydt(15) = a_APCm_TNFa * (1 + alpha_IL1b_TNFa * IL1b_PT / ( IL1b_PT + K_IL1b_TNFa)) * K_IL10_TNFa ^ hill_IL10_TNFa / ( IL10_PT ^ hill_IL10_TNFa + K_IL10_TNFa ^ hill_IL10_TNFa)...
%     * K_CRY ^ 1 / (K_CRY ^ 1 + CRY ^ 1) * APCm_PT^hill_APCm_TNFa/(APCm_PT^hill_APCm_TNFa + K_APCm_TNFa^hill_APCm_TNFa) + ...
%     a_Th1_TNFa * Th1_PT + a_CTL_TNFa * CTL_PT + gamma_BPT * c_BPT_TNFa * TNFa_B - c_PTLN_TNFa * TNFa_PT - c_PTB_cyto * TNFa_PT - d_TNFa * TNFa_PT;

% TNFa
if control_mphi == 1     %mphi使用默认的稳态生物钟
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
% dydt(15) = a_APCm_TNFa*2.6*10^8*APCm_PT^5 + ...
%     a_Th1_TNFa * Th1_PT + a_CTL_TNFa * CTL_PT + gamma_BPT * c_BPT_TNFa * TNFa_B - c_PTLN_TNFa * TNFa_PT - c_PTB_cyto * TNFa_PT - 10 * TNFa_PT;

% IFNg
dydt(16) = a_Neu_IFNg * Neu_PT + a_CTL_IFNg * CTL_PT + gamma_BPT * c_BPT_cyto * IFNg_B - c_PTLN_cyto * IFNg_PT - c_PTB_cyto * IFNg_PT - d_IFNg * IFNg_PT;

% IL1b不需要
dydt(17) = 0;

% IL2
dydt(18) = a_Th1_IL2 * Th_PT + gamma_BPT * c_BPT_cyto * IL2_B - c_PTLN_cyto * IL2_PT - c_PTB_cyto * IL2_PT - d_IL2 * IL2_PT;

% IL4
dydt(19) = a_Th2_IL4 * Th_PT + gamma_BPT * c_BPT_cyto * IL4_B - c_PTLN_cyto * IL4_PT - c_PTB_cyto * IL4_PT - d_IL4 * IL4_PT;

% IL6
if control_mphi == 1     %mphi使用默认的稳态生物钟
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

% dydt(20) = a_APCm_IL6 *100000000* APCm_PT^4 + a_Th17_IL6 * Th17_PT + gamma_BPT * c_BPT_IL6 * IL6_B - c_PTLN_IL6 * IL6_PT - c_PTB_cyto * IL6_PT - d_IL6 * IL6_PT;

% IL10
if control_mphi == 1     %mphi使用默认的稳态生物钟
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
alpha_CCL1921 = CCL21_LN ^ 2 / (K_CCL21 ^ 2 + CCL21_LN ^ 2) + alpha_CCL19 * CCL19_LN ^ 2 / (K_CCL19 ^ 2 + CCL19_LN ^ 2);
alpha_CCL1921_DC = CCL21_LN ^ 2 / (K_CCL21_DC ^ 2 + CCL21_LN ^ 2) + alpha_CCL19_DC * CCL19_LN ^ 2 / (K_CCL19_DC ^ 2 + CCL19_LN ^ 2);
alpha_egress = S1PR1 ^ n_S1PR1 / (S1PR1 ^ n_S1PR1 + K_S1PR1 ^ n_S1PR1);

% DC_LN
% dydt(23) = gamma_PTLN * b * c_PTLN_cell * CCL1921_DC_LN ^ n_CCL1921_APCm / (K_CCL1921_APCm ^ n_CCL1921_APCm + CCL1921_DC_LN ^ n_CCL1921_APCm) * DC_PT - d_DC * DC_LN;
dydt(23) = gamma_PTLN * c_PTLN_cell * alpha_CCL1921_DC * DC_PT - d_DC * DC_LN;
% gamma_PTLN * c_PTLN_cell * DC_PT * ( CCL21 ^ 2 / (K_CCL21_DC ^ 2 + CCL21 ^ 2) + alpha * CCL19 ^ n_C / (K_CCL19_DC ^ n_C + CCL19 ^ n_C)) - d_DC * DC_LN;

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
% dydt(26) = a_CD4Ta_Th *(1 + alpha_IL2_Th * IL2_LN / (IL2_LN + K_IL2_Th) + alpha_IL6_Th * IL6_LN / (IL6_LN + K_IL6_Th) + alpha_IFNg_Th * IFNg_LN / (IFNg_LN + K_IFNg_Th)) * ...
%     K_IL10_Th / (K_IL10_Th + IL10_LN) * K ^ 2 / (K ^ 2 + DC_LN ^ 2) * CD4Ta_LN - d_Treg_T * Th_LN * Treg_LN + gamma_BLN * alpha_HEC_CD4T * c_BLN_cell * alpha_CCL1921 * Th_B...
%     - c_LNB_cell * alpha_egress * Th_LN - d_Th * Th_LN + gamma_PTLN * c_PTLN_cell * Th_PT;
dydt(26) = a_CD4Ta_Th * 2 * REV / (REV + K_REV_Th) * (1 + alpha_IL2_Th * IL2_LN / (IL2_LN + K_IL2_Th) + alpha_IL6_Th * IL6_LN / (IL6_LN + K_IL6_Th) + alpha_IFNg_Th * IFNg_LN / (IFNg_LN + K_IFNg_Th)) * ...
    K_IL10_Th / (K_IL10_Th + IL10_LN) * K ^ 2 / (K ^ 2 + DC_LN ^ 2) * CD4Ta_LN - d_Treg_T * Th_LN * Treg_LN + gamma_BLN * alpha_HEC_CD4T * c_BLN_cell * alpha_CCL1921 * Th_B...
    - c_LNB_cell * alpha_egress * Th_LN - d_Th * Th_LN + gamma_PTLN * c_PTLN_cell * Th_PT;

% dydt(26) = a_CD4Ta_Th1 * (1 + alpha_IFNg_Th1 * IFNg_LN / (IFNg_LN + K_IFNg_Th1)) * K_IL4_Th1 / (K_IL4_Th1 + IL4_LN) * K_IL10_Th1 / (K_IL10_Th1 + IL10_LN) * CD4Ta_LN + ...
%     gamma_BLN * c_BLN_Th1 * (1 + alpha_TNFa_Th1 * TNFa_LN / (K_TNFa_Th1 + TNFa_LN)) * (CCL1921_LN ^ 2 / (K_CCL1921_Th1 ^ 2 + CCL1921_LN ^ 2)) * Th1_B...
%     - c_LNB_Th1 * S1PR1 ^ 2 / (S1PR1 ^ 2 + K_S1PR1_Th1 ^ 2) * Th1_LN - d_Th * Th1_LN - d_Treg_T * Th1_LN * Treg_LN + gamma_PTLN * c_PTLN_Th1 * Th1_PT;


% hHEC(proliferation via VEGF signal pathway)
% dydt(27) = -c_hHEC_iHEC * (DC_LN ^ n_DC_HEC / (DC_LN ^ n_DC_HEC + K_APCm_iHEC ^ n_DC_HEC) + alpha_iHEC * TNFa_LN / (TNFa_LN + K_TNFa_iHEC)) * hHEC + ...
%     c_iHEC_hHEC * iHEC + rho_HEC * hHEC * (1 - (hHEC + iHEC) / N_HEC) * (1 + alpha_APCm_iHEC * DC_LN ^ n_DC_HEC / (DC_LN ^ n_DC_HEC + K_APCm_iHEC ^ n_DC_HEC))...
%     - d_HEC * hHEC;
dydt(27) = -c_hHEC_iHEC * (DC_LN ^ n_DC_HEC / (DC_LN ^ n_DC_HEC + K_APCm_iHEC ^ n_DC_HEC) + alpha_iHEC * TNFa_LN / (TNFa_LN + K_TNFa_iHEC)) * hHEC + ...
    c_iHEC_hHEC * iHEC + rho_HEC * hHEC * (1 - (hHEC + iHEC) / N_HEC) ...
    - d_HEC * hHEC;
%iHEC(inflamed-state transition via VEGF and TNF receptor)
% dydt(28) = c_hHEC_iHEC * (DC_LN ^ n_DC_HEC / (DC_LN ^ n_DC_HEC + K_APCm_iHEC ^ n_DC_HEC) + alpha_iHEC * TNFa_LN / (TNFa_LN + K_TNFa_iHEC)) * hHEC + ...
%     rho_HEC * iHEC * (1 - (hHEC + iHEC) / N_HEC) * (1 + alpha_APCm_iHEC * DC_LN ^ n_DC_HEC / (DC_LN ^ n_DC_HEC + K_APCm_iHEC ^ n_DC_HEC))...
%     - c_iHEC_hHEC * iHEC - d_HEC * iHEC;
dydt(28) = c_hHEC_iHEC * (DC_LN ^ n_DC_HEC / (DC_LN ^ n_DC_HEC + K_APCm_iHEC ^ n_DC_HEC) + alpha_iHEC * TNFa_LN / (TNFa_LN + K_TNFa_iHEC)) * hHEC + ...
    rho_HEC * iHEC * (1 - (hHEC + iHEC) / N_HEC) ...
    - c_iHEC_hHEC * iHEC - d_HEC * iHEC;
% CCL1921_DC
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
% dydt(37) = DC_LN / (DC_LN + K_APCm_BGC) * (1 ) * ...
%     (a_Bn_BGC * alpha_activation_CB * Bn_LN + a_Bm_BGC * Bm_LN + rho_BGC * DI * alpha_Tfh_BGC * 0.5 * BGC_LN * (1 - BGC_LN / N_GC)) ...
%     - a_BGC_Bm * BGC_LN - d_BGC * BGC_LN;
% ASC

dydt(38) = gamma_BLN * c_BLN_cell * alpha_HEC * alpha_CCL1921 * ASC_B...
    - c_LNB_cell * alpha_egress * ASC_LN + gamma_PTLN * c_PTLN_cell * ASC_PT + a_BGC_ASC * BGC_LN * K ^ 2 / (K ^ 2 + DC_LN ^ 2)...
    - d_ASC * ASC_LN;

% Bm
dydt(39) = a_BGC_Bm * BGC_LN - a_Bm_BGC * (1 + alpha_IL4_BGC * IL4_LN / (K_IL4_BGC + IL4_LN) + alpha_IL6_BGC * IL6_LN / (K_IL6_BGC + IL6_LN)) * ...
    DC_LN / (DC_LN + K_APCm_BGC) * Bm_LN - d_Bm * Bm_LN;

% CCL21
if control_Lymphocyte_rhythmic_homing == 1     %EC使用默认的稳态生物钟
    dydt(40) = a_0_CCL1921 * (CLOCK_BMAL1 .^ 3 ./ (K_BC_CCL1921 .^ 3 + CLOCK_BMAL1 .^ 3) ) - d_CCL1921 * CCL21_LN;
else
    dydt(40) = a_0_CCL1921 * 0.467 - d_CCL1921 * CCL21_LN; % 否则EC使用动态行为下的特殊生物钟（敲除，0态出发，无振幅等）
end

% S1PR1
if control_Lymphocyte_rhythmic_homing == 1     %T使用默认的稳态生物钟
    dydt(41) = a_S1PR1 * (K_BC_S1PR1 .^ 3 ./ (K_BC_S1PR1 .^ 3 + (CLOCK_BMAL1_t .^ 3) )) - d_S1PR1 * S1PR1;
else
    dydt(41) = a_S1PR1 * 0.0634 - d_S1PR1 * S1PR1; % 否则T使用动态行为下的特殊生物钟（敲除，0态出发，无振幅等）
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
% CCL19
dydt(50) = a * a_APCm_CCL1921 * DC_LN ^ n_sec_CCL19 / (DC_LN ^ n_sec_CCL19 + K_DC_CCL1921 ^ n_sec_CCL19) - d_CCL1921 * CCL19_LN;

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
function dydt = ode_DCCCL_time_reversal(t, y, Par, Time_circadian_clock, ZT_initial)
circa_t_index = mod(floor((t * 24 + ZT_initial) / 24 * 1000), 1000) + 1; 
CLOCK_BMAL1_t = Time_circadian_clock(circa_t_index, 7); 

CCL19 = y(1); 
DC_LN = y(2); 
DC_PT = y(3); 
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

dydt = zeros(3, 1); 

% soluble CCL21 
dydt(1) = a_DC_CCL * DC_LN ^ n_D / (DC_LN ^ n_D + K_DC_CCL ^ n_D) - d_CCL * CCL19; % x-CCL
% DC in LN
dydt(2) = gamma_PTLN * c_PTLN_cell * DC_PT * (alpha * CCL19 ^ n_C / (K_CCL19_DC ^ n_C + CCL19 ^ n_C)) - d_DC * DC_LN; %y-DC
dydt(3) = 0; 
% time_reverse
dydt = -dydt; 
end
function dydt = ode_DCCCL_no_time_reverse(t, y, Par, Time_circadian_clock, ZT_initial)
circa_t_index = mod(floor((t * 24 + ZT_initial) / 24 * 1000), 1000) + 1; 
CLOCK_BMAL1_t = Time_circadian_clock(circa_t_index, 7); 

CCL19 = y(1); 
DC_LN = y(2); 
DC_PT = y(3); 
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

dydt = zeros(3, 1); 

% soluble CCL21 and CCL19
dydt(1) = a_DC_CCL * DC_LN ^ n_D / (DC_LN ^ n_D + K_DC_CCL ^ n_D) - d_CCL * CCL19; % x-CCL
% DC in LN
dydt(2) = gamma_PTLN * c_PTLN_cell * DC_PT * (alpha * CCL19 ^ n_C / (K_CCL19_DC ^ n_C + CCL19 ^ n_C)) - d_DC * DC_LN; %y-DC
dydt(3) = 0; 
% no time_reverse
% dydt=-dydt;
end
