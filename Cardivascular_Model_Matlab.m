clc ; clear all ; close all ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Simulation of the Cardiovascular System %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Flags:
Plot_flag = 1; % 0 = off , 1 = on
 
%% system parameters:
% Time parameters
HR           = 60 + 9+7+9 ; % [BPM] % 60 + sum of last digits from all members
dt           = 5e-4        ; % [sec]
Heart_cycles = 20          ; % total heart cycles
N_per_cycle      = round(60/(HR*dt))    ; % number of steps per heart cycle
  
% Heart parameters
V0           = 15  ; % [ml] V0 for Plv calculation
Emax         = 2.0 ; % contractility
N_Systole    = round((1/3) * N_per_cycle) ; % number of points per ventricle Systole
N_Diastole   = round((2/3) * N_per_cycle); % number of points per ventricle Diastole
En(1:N_Systole) = 0.5*(1 + sin((2*pi*(1:N_Systole)/N_Systole)-(pi/2)));
En(N_Systole+1:N_per_cycle) = 0;
E            = max((10/105), Emax*En); % Heart Elasticity: combine Systole (elasticity function) & Diastole (uniform value) Ed=10/105; %diastolic E (constant value) 
 
% Vascular constants:
Ra = 0.1;  % arterial resistance 
Rp = 1.0;  % peripheral resistance
Rv = 0.01; % venous filling resistance
Ca = 2.0;  % arterial compliance
Cv = 300.0; % venous compliance 
 
% Initiate variables:
%Volume [ml]
Vlv(1)  = 120;  % left ventricle
Va(1)   = 270;  % arteries
Vv(1)   = 2700; % veins 
BV_100=Vlv(1)+Va(1)+Vv(1);
BV_axis=linspace(BV_100*0.5,BV_100*3); %Blood Volume for Q4
BV_axis_norm=BV_axis./BV_100;
%Pressure [mmHg]
Plv(1)  = 0;    % left ventricle
Pa(1)   = 70;   % arterial capacitor
Pv(1)   = 9;    % venous filling 
Pao(1)  = 100;   % aorta
%Flow [ml/sec]
Qlv(1)  = 0;    % left ventricle (outflow)
Qp(1)   = 0;    % peripheral resistance
Qv(1)   = 0;    % ventricle filling (inflow)
avg_Pao=[];
%% Main Program
for CycleIdx = 1 : Heart_cycles % main loop for each heart cycle
    for StepInCycle = 2 : N_per_cycle 
	%calculating all variables for each cycle at N points:
        %Volumes [ml]
		Vlv(StepInCycle) = Vlv(StepInCycle-1) + (Qv(StepInCycle-1) - Qlv(StepInCycle-1))*dt;
        Va(StepInCycle) = Va(StepInCycle-1) + (Qlv(StepInCycle-1) - Qp(StepInCycle-1))*dt;
        Vv(StepInCycle) = Vv(StepInCycle-1) + (Qp(StepInCycle-1) - Qv(StepInCycle-1))*dt;
        %Pressures [mmHg] ...
        Plv(StepInCycle) = E(StepInCycle)*(Vlv(StepInCycle)-V0);
        Pa(StepInCycle) = Va(StepInCycle)/Ca;
        Pv(StepInCycle) = Vv(StepInCycle)/Cv;
        %Flows [ml/sec] ...
        Qv(StepInCycle) = max(0,((Pv(StepInCycle)-Plv(StepInCycle))/Rv));
        Qlv(StepInCycle) = max(0,((Plv(StepInCycle)-Pa(StepInCycle))/Ra));
        Qp(StepInCycle) = ((Pa(StepInCycle) - Pv(StepInCycle))/Rp);
        
        Pao(StepInCycle) = Pa(StepInCycle) + Ra*Qlv(StepInCycle);
        
    end

    % Min/Max of Pao
    % MAX
    [M,I] = max(Pao); Pao_max(CycleIdx) = M;
    Pao_max_index(CycleIdx) = dt*(I + ((CycleIdx-1)*N_per_cycle));
    % MIN
    [M,I] = min(Pao); Pao_min(CycleIdx) = M;
    Pao_min_index(CycleIdx) = dt*(I + ((CycleIdx-1)*N_per_cycle));

    % Save each variable from the current cycle to a continuos variable:
    %VolumeS [ml] 

    c_Vlv(CycleIdx*length(Vlv)-(length(Vlv)-1):CycleIdx*length(Vlv)) = Vlv;
    c_Va(CycleIdx*length(Va)-(length(Va)-1):CycleIdx*length(Va)) = Va;
    c_Vv(CycleIdx*length(Vv)-(length(Vv)-1):CycleIdx*length(Vv)) = Vv;
    
    %PressureS [mmHg] 

    c_Plv(CycleIdx*length(Plv)-(length(Plv)-1):CycleIdx*length(Plv)) = Plv;
    c_Pa(CycleIdx*length(Pa)-(length(Pa)-1):CycleIdx*length(Pa)) = Pa;
    c_Pao(CycleIdx*length(Pao)-(length(Pao)-1):CycleIdx*length(Pao)) = Pao;
    c_Pv(CycleIdx*length(Pv)-(length(Pv)-1):CycleIdx*length(Pv)) = Pv;
    
    %FlowS [ml/sec] 
    
    c_Qlv(CycleIdx*length(Qlv)-(length(Qlv)-1):CycleIdx*length(Qlv)) = Qlv;
    c_Qp(CycleIdx*length(Qp)-(length(Qp)-1):CycleIdx*length(Qp)) = Qp;
    c_Qv(CycleIdx*length(Qv)-(length(Qv)-1):CycleIdx*length(Qv)) = Qv;

    % Update the initial variables before the next cycle: 
    %Volume [ml] 
    Vlv(1) = Vlv(end); Va(1) = Va(end); Vv(1) = Vv(end);
    %Pressure [mmHg] 
    Plv(1) = Plv(end); Pa(1) = Pa(end); Pao(1) = Pao(end); Pv(1) = Pv(end);
    %Flow [ml/sec] 
	Qlv(1) = Qlv(end); Qp(1) = Qp(end); Qv(1) = Qv(end);
    
end

%% Plots

if Plot_flag
    %1.1
    figure(1); %Aortic Pressure
    t = (1:length(c_Pao))*dt; % Time vector
    plot(t, c_Pao); hold on;         % Graph 1: Aortic Pressure
    scatter(Pao_max_index,Pao_max);  % Scatter max pressure points
    scatter(Pao_min_index,Pao_min);  % Scatter min pressure points
    title('Aortic Pressure');
    xlabel('Time [sec]'); ylabel('Pressure [mmHg]');
    hold off
    
    %2
    figure(2); 
    subplot(2,1,1) %One Heart cycle
    min_indices=Pao_min_index/dt; %Calculate one cycle between minimum pressure points 
    sliced_idx=min_indices(17)-400:min_indices(18)+400; %slicing one heart cycle
    sliced_t=sliced_idx*dt; % time vector
    plot(sliced_t,c_Pa(sliced_idx),'MarkerSize',8); hold on;
    plot(sliced_t,c_Pao(sliced_idx),'MarkerSize',8);
    plot(sliced_t,c_Plv(sliced_idx),'MarkerSize',8);
    plot(sliced_t,c_Pv(sliced_idx),'MarkerSize',8);
    grid on;
    title('One Heart Cycle: Pressure vs. Time')
    xlabel('Time [Sec]')
    ylabel('Pressure [mmHg]')
    legend('Arterial Pressure','Aortic Pressure','Left Ventricle Pressure','Veins Pressure');
    text([11.3115, 11.3415, 11.4515, 11.51 ],[9.4881, 64.4917, 83.1983, 9.53133],{'\rightarrow Left AV Valve Closes','\rightarrow Aortic Valve Opens',  '\leftarrow  Aortic Valve Closes', '\leftarrow  Left AV valve Opens'});
    
    subplot(2,1,2) %Left Ventricle volume
    plot(sliced_t,c_Vlv(sliced_idx)) ; hold on
    xlabel('Time [Sec]')
    ylabel('Volume [ml]')
    title('Left Ventricle Volume')
    [EDV,EDV_idx]=max(c_Vlv(sliced_idx)); %EDV
    [ESV,ESV_idx]=min(c_Vlv(sliced_idx)); %ESV
    text((EDV_idx*dt)+((min_indices(17)-400)*dt),EDV,'End Diastolic Volume')
    text((ESV_idx*dt)+((min_indices(17)-400)*dt),ESV,'End systolic Volume')    
    
    %Plot SV annotation
    SV = [ESV EDV]; 
    SV_idx = (EDV_idx*dt)+((min_indices(17)-400)*dt);  
    plot([SV_idx SV_idx],SV,'--k','Linewidth',1.4)
    text(SV_idx,(ESV+EDV)/2,'\leftarrow Stroke Volume')
    grid on;
    
    %% Q3.2
    figure(3);
    %t2={114.115;9.4881;0}; t3={114.115;64.4917;0.03} ; t4={69.2253;83.1983;0.14} ; t1={69.2253;9.5313;0.1985};
    Volume_t=[114.115 114.115 69.2253 69.2253];
    Pressure_t=[9.4881 64.4917 83.1983 9.5313];
    time_e=0.0165+[0 0.03 0.14 0.1985];
    elasticity=Pressure_t./Volume_t;               %Calculate elasticity based on 4 points
    m1=(elasticity(2)-elasticity(1))/(0.03-0);
    m2=(elasticity(3)-elasticity(2))/(0.14-0.03);   %Defining the middle line
    m3=(elasticity(4)-elasticity(3))/(0.1985-0.14);
    b1=elasticity(2)-0.03*m1;
    b2=elasticity(2)-0.03*m2;                   %using linear equation
    b3=elasticity(3)-0.14*m3;
    new_t=0.0165+[linspace(0,0.03,1000)  linspace(0.03 ,0.14,1000) linspace(0.14,0.1985,1000)];
    new_e=[m1*linspace(0,0.03,1000)+b1  m2*linspace(0.03 ,0.14,1000)+b2 m3*linspace(0.14,0.1985,1000)+b3];
    plot(new_t, new_e) ; hold on 
    sin_E_time = (0:N_per_cycle-1)*dt; %time vector
    plot(sin_E_time,E)
    plot((new_t(1000)+new_t(2000))/2,(new_e(1000)+new_e(2000))/2 ,'o-k', 'MarkerFaceColor', 'r','MarkerSize',8)
    title('Heart Elasticity over Time');
    xlabel('Time(single heart cycle) [sec]'); ylabel('Elasticity [mmHg/mL]');
    legend('P-V Loop Elastance','sin Elastance','t_d'); 
    
%% Q4
for i=1:length(BV_axis) %Initate simulation for BV testing
 % Initiate variables:
    %Volume [ml]
    Vlv(1)  = 120*BV_axis_norm(i);  % left ventricle
    Va(1)   = 270*BV_axis_norm(i);  % arteries
    Vv(1)   = 2700*BV_axis_norm(i); % veins 
    %Pressure [mmHg]
    Plv(1)  = 0;    % left ventricle
    Pa(1)   = 70;   % arterial capacitor
    Pv(1)   = 9;    % venous filling 
    Pao(1)  = 100;   % aorta
    %Flow [ml/sec]
    Qlv(1)  = 0;    % left ventricle (outflow)
    Qp(1)   = 0;    % peripheral resistance
    Qv(1)   = 0;    % ventricle filling (inflow)
   
    %% Main Program
    for CycleIdx = 1 : Heart_cycles % main loop for each heart cycle
        for StepInCycle = 2 : N_per_cycle 
            %Volumes [ml]
		    Vlv(StepInCycle) = Vlv(StepInCycle-1) + (Qv(StepInCycle-1) - Qlv(StepInCycle-1))*dt;
            Va(StepInCycle) = Va(StepInCycle-1) + (Qlv(StepInCycle-1) - Qp(StepInCycle-1))*dt;
            Vv(StepInCycle) = Vv(StepInCycle-1) + (Qp(StepInCycle-1) - Qv(StepInCycle-1))*dt;
            %Pressures [mmHg] ...
            Plv(StepInCycle) = E(StepInCycle)*(Vlv(StepInCycle)-V0);
            Pa(StepInCycle) = Va(StepInCycle)/Ca;
            Pv(StepInCycle) = Vv(StepInCycle)/Cv;
            %Flows [ml/sec] ...
            Qv(StepInCycle) = max(0,((Pv(StepInCycle)-Plv(StepInCycle))/Rv));
            Qlv(StepInCycle) = max(0,((Plv(StepInCycle)-Pa(StepInCycle))/Ra));
            Qp(StepInCycle) = ((Pa(StepInCycle) - Pv(StepInCycle))/Rp);
            Pao(StepInCycle) = Pa(StepInCycle) + Ra*Qlv(StepInCycle);
            
        end
        %VolumeS [ml] ...
        c_Vlv(CycleIdx*length(Vlv)-(length(Vlv)-1):CycleIdx*length(Vlv)) = Vlv;
        c_Va(CycleIdx*length(Va)-(length(Va)-1):CycleIdx*length(Va)) = Va;
        c_Vv(CycleIdx*length(Vv)-(length(Vv)-1):CycleIdx*length(Vv)) = Vv;
        %PressureS [mmHg] 
        c_Plv(CycleIdx*length(Plv)-(length(Plv)-1):CycleIdx*length(Plv)) = Plv;
        c_Pa(CycleIdx*length(Pa)-(length(Pa)-1):CycleIdx*length(Pa)) = Pa;
        c_Pao(CycleIdx*length(Pao)-(length(Pao)-1):CycleIdx*length(Pao)) = Pao;
        c_Pv(CycleIdx*length(Pv)-(length(Pv)-1):CycleIdx*length(Pv)) = Pv;
        %FlowS [ml/sec] 
        c_Qlv(CycleIdx*length(Qlv)-(length(Qlv)-1):CycleIdx*length(Qlv)) = Qlv;
        c_Qp(CycleIdx*length(Qp)-(length(Qp)-1):CycleIdx*length(Qp)) = Qp;
        c_Qv(CycleIdx*length(Qv)-(length(Qv)-1):CycleIdx*length(Qv)) = Qv;
        %Volume [ml] ...
        Vlv(1) = Vlv(end); Va(1) = Va(end); Vv(1) = Vv(end);
        %Pressure [mmHg] ...
        Plv(1) = Plv(end); Pa(1) = Pa(end); Pao(1) = Pao(end); Pv(1) = Pv(end);
        %Flow [ml/sec] ...
	    Qlv(1) = Qlv(end); Qp(1) = Qp(end); Qv(1) = Qv(end);
    end
    avg_Pao(end+1) =mean(c_Pao(14*N_per_cycle:end)); %average Pao while stable 
end
figure(4); %Plotting results
plot(BV_axis,avg_Pao); hold on;
plot(BV_axis(21),avg_Pao(21),'o-k', 'MarkerFaceColor', 'r','MarkerSize',8);
title('Average stable Aortic Pressure over Blood Volume');
xlabel('Blood Volume [ml]');
ylabel('Aortic Pressure[mmHg]');
xlim([BV_axis(1),BV_axis(end)]);
legend('Pressure vs BV','100% BV marker')
end


