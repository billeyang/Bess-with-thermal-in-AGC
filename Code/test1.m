%huodian_data = xlsread("火储典型.xlsx");
clc;
clear all;
AGC_ask = xlsread("火储典型.xlsx",'N2:N601'); %读取AGC曲线
%POWER_REF = xlsread("火储典型2-t.xlsx",'G2:G601');%读取火电出力曲线
POWER_OUT = xlsread("火储典型.xlsx",'O2:O601');%读取火电前移t时刻出力曲线，主要是因为火电实际反应要比AGC指令慢，因此按照t时刻前移对应火电出力
x = AGC_ask - POWER_OUT; %计算AGC和火电出力差值
ess = zeros(600,1); %建立储能荷电状态参数，起始值为0
%total_charge = zeros(600,1);
%total_discharge = zeros(600,1);
 
power_constraint = 36;%储能功率36MW
energy_capacity = 18; %储能能量18
charge_ratio = 0.9;%充电深度0.9
discharge_ratio = 0.1;%放电深度0.1
%energy_constraint = 7.2;
energy_constraint = energy_capacity *(charge_ratio-discharge_ratio);%实际储能能量
time_period = 0.04; % 2min24s/60min = 144s/3600s = 0.04


%ess(i)为储能当前能量
for i = 1 : 600 %该电厂一天调用次数为600
    %a(i) = Diff_ori(i);
    %% 放电-AGC与出力差值小于储能功率
    if x(i) < power_constraint & x(i) > 0 
        if ess(i) >  x(i) * time_period %当电池储存的能量大于差值需要的能量时，放出差值*时间的能量
            ess(i+1) = ess(i) - x(i) * time_period;
            %total_discharge(i+1) = total_discharge(i) + x(i) * time_period;
        elseif ess(i) >0 && ess(i)< (x(i) * time_period)%当电池储存的能量小于差值需要的能量时，放出电池的全部能量
            ess(i+1) = 0;
            %total_discharge(i+1) = total_discharge(i) + ess(i);
        elseif ess(i) == 0%当电池储存的能量为0，无动作。
            ess(i+1) = ess(i);
        end
    %% 放电-AGC与出力差值大于储能功率
    elseif x(i) >= power_constraint
        if ess(i) >  power_constraint * time_period%当电池储存的能量大于PCS功率*时间需要的能量时，放出功率*t时间的能量
            ess(i+1) = ess(i) - power_constraint * time_period;
            %total_discharge(i+1) = total_discharge(i) + power_constraint * time_period;
        elseif ess(i) >0 && ess(i)< (power_constraint * time_period)%当电池储存的能量小于PCS功率*t时间需要的能量时，放出电池的全部能量
            ess(i+1) = 0;
            %total_discharge(i+1) = total_discharge(i) + ess(i);
        elseif ess(i) == 0%当电池储存的能量为0，无动作。
            ess(i+1) = ess(i);
            %total_discharge(i+1) = total_discharge(i);
        end
        %% 充电-AGC与出力差值小于储能功率的绝对值
    elseif x(i) < 0 & x(i) > -power_constraint
        if (ess(i) - x(i)* time_period) <  energy_constraint%当电池储存的能量+下一个时段充电的能量小于能量约束，充入差值*t时间的能量
            ess(i+1) = ess(i) - x(i) * time_period;
            %total_charge(i+1) = total_charge(i) - x(i) * time_period;
        elseif ess(i) <= energy_constraint && (ess(i) + (power_constraint * time_period)) > energy_constraint%当电池储存的能量小于能量约束，同时存储的能量+下一个时段充电的能量大于能量约束，充满
            %elseif (ess(i) - (0.04*x(i))) > energy_constraint
            ess(i+1) = energy_constraint;
            %total_charge(i+1) = total_charge(i) + energy_constraint - ess(i);
        elseif ess(i) == energy_constraint %当电池充满，无动作。
            ess(i+1) = ess(i);
            %total_charge(i+1) = total_charge(i);
        end
        
        %% 充电-AGC与出力差值大于储能功率的绝对值        
    elseif x(i) <= -power_constraint
        if (ess(i) + power_constraint * time_period) <  energy_constraint%当电池储存的能量+下一个时段功率约束*t小于能量约束，充入功率约束*t的能量
            ess(i+1) = ess(i) + power_constraint * time_period;
            %total_charge(i+1) = total_charge(i) + power_constraint * time_period;
        elseif ess(i) <= energy_constraint && (ess(i) + (power_constraint * time_period)) > energy_constraint%当电池储存的能量小于能量约束，同时存储的能量+功率约束*t的能量大于能量约束，充满
            %elseif (ess(i) - (0.04*x(i))) > energy_constraint
            ess(i+1) = energy_constraint;
            %total_charge(i+1) = total_charge(i) + energy_constraint - ess(i);
        elseif ess(i) == energy_constraint%当电池充满，无动作。
            ess(i+1) = ess(i);
            %total_charge(i+1) = total_charge(i);
        end
        %% 没电-当差值=0，无动作
    elseif x(i) ==0 
        ess(i+1) = ess(i);
        %total_charge(i+1) = total_charge(i);
    end
    
    bess_act(i+1,1) = ess(i+1) - ess(i); %每个时段储能的能量
    SOC(i,1) = 0.1+(ess(i)/energy_capacity);%每个时段储能的荷电状态
    bess_power_act(i + 1,1) = bess_act(i + 1)/0.04; %每个时段储能提供的功率支持
    
end
for i = 1: 601
    if bess_act(i) >= 0
        bess_charge(i,1) = bess_act(i);%每个时段储能的充电电量
        bess_discharge(i,1) = 0;
    else
        bess_charge(i,1) = 0;
        bess_discharge(i,1) = bess_act(i);%每个时段储能的放电电量
    end
    
end
BESS_OUT =  POWER_OUT(1:600) - bess_power_act(2:601);

%% excel导出
xlswrite('output.xlsx', ess(2:601), 1, 'C2');%输出储能的当前能量
xlswrite('output.xlsx', x, 1, 'D2');%输出火电AGC和火电本体的功率差值
xlswrite('output.xlsx', AGC_ask, 1, 'E2');%输出AGC原始需求
xlswrite('output.xlsx', POWER_OUT, 1, 'F2');%输出火电出力需求
xlswrite('output.xlsx', BESS_OUT, 1, 'G2');%输出火电+储能出力
xlswrite('output.xlsx', bess_power_act(2:601), 1, 'H2');%输出储能时段功率
xlswrite('output.xlsx', bess_charge, 1, 'I2');%输出储能时段充电能量
xlswrite('output.xlsx', bess_discharge, 1, 'J2');%输出储能时段放电能量

figure;
title('AGC/纯火电/火储联调三种出力模式对比')
xlabel('时间0~24h(每个点间距2min24s)')
ylabel('功率(单位：MW)')
hold on;
plot(AGC_ask,'g');
grid on;
hold on;
plot(POWER_OUT,'r:');
hold on;
% plot(POWER_REF,'c-.');
% hold on;
plot(BESS_OUT,'b--');
legend('AGC指令','纯火电','火储联调','location','northwest');
figure;
plot(SOC,'--r','MarkerEdgeColor','b');
title('储能系统SOC情况')
xlabel('时间0~24h(每个点间距2min24s)')
ylabel('SOC')
legend('SOC曲线','location','northwest');
figure;
plot(bess_power_act(2:601),'-.g','MarkerEdgeColor','b');
title('储能本体实时出力')
xlabel('时间0~24h(每个点间距2min24s)');
ylabel('储能功率');
legend('储能本体出力曲线','location','northwest');
