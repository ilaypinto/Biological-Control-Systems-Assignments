function output = Pao_PID(Kp,Kd,Ki)
%% Setting default PID parameters as 0
if nargin==0
    Kp=0; Kd=0; Ki=0;
elseif nargin==1
    Kd=0; Ki=0;
elseif nargin==2
    Ki=0;
end
%% system parameters: Taken directly from project's part A
% Time parameters
HR           = 80 ; % [BPM] 80 as default, as stated in the assignment
dt           = 5e-4        ; % [sec]
Heart_cycles = 400          ; % total heart cycles
N_per_cycle      = round(60/(HR*dt))    ; % number of steps per heart cycle
  
% Heart parameters
V0           = 15  ; % [ml] V0 for Plv calculation
Emax         = 2.0 ; % contractility
N_Systole    = round((1/3) * N_per_cycle) ; % number of points per ventricle Systole
N_Diastole   = round((2/3) * N_per_cycle); % number of points per ventricle Diastole
En(1:N_Systole) = 0.5*(1 + sin((2*pi*(1:N_Systole)/N_Systole)-(pi/2)));
En(N_Systole+1:N_per_cycle) = 0;
E = max((10/105), Emax*En); % Heart Elasticity: combine Systole (elasticity function) & Diastole (uniform value) Ed=10/105; %diastolic E (constant value) 
 
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

%Pressure [mmHg]
Plv(1)  = 0;    % left ventricle
Pa(1)   = 70;   % arterial capacitor
Pv(1)   = 9;    % venous filling 
Pao(1)  = 100;  % aorta

%Flow [ml/sec]
Qlv(1)  = 0;    % left ventricle (outflow)
Qp(1)   = 0;    % peripheral resistance
Qv(1)   = 0;    % ventricle filling (inflow)

%Pao_in=79.8677745032775; %Taken from Avg Pao calc in part A of the project for HR 85
%Pao_in=78.2237923022584;%For HR 80
Pao_in=77.828;
error = zeros(1,Heart_cycles);
error_sum=0;
output=[];

%% Main Program
for CycleIdx = 1 : Heart_cycles % main loop for each heart cycle
    % Bleeding 15% of Blood Volume for cycle 70
    if (CycleIdx == 200)
        Vlv(1)=Vlv(1)*0.85; 
        Va(1)=Va(1)*0.85;
        Vv(1)=Vv(1)*0.85;
    end
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

    Pao_out=mean(Pao(1:N_per_cycle));
    output=[output Pao_out];

    % Update the initial variables before the next cycle: 
    %Volume [ml] 
    Vlv(1) = Vlv(end); Va(1) = Va(end); Vv(1) = Vv(end);
    %Pressure [mmHg] 
    Plv(1) = Plv(end); Pa(1) = Pa(end); Pao(1) = Pao(end); Pv(1) = Pv(end);
    %Flow [ml/sec] 
	Qlv(1) = Qlv(end); Qp(1) = Qp(end); Qv(1) = Qv(end);

    %error
    error(CycleIdx) = Pao_in-Pao_out;
    error_sum = error_sum+error(CycleIdx);
    %calculation the PID
    P = error(CycleIdx)*Kp;
    I = error_sum*Ki;
    if (CycleIdx ~= 1)
        D = (error(CycleIdx)-error(CycleIdx-1))*Kd;
    else
        D = error(CycleIdx)*Kd;
    end
    PID = P + I + D + 1;

    %New Parameters using PID:
    Emax = max(0,Emax*PID);
    E = max((10/105),Emax*En);
    Cv = max(0.001,Cv*PID);
    Rp = max(0.001,Rp*PID);
    
end
end