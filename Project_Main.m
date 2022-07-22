%% Main code for Matlab project part B
%Run code with Pao_PID.m in the same file. Don't run Pao_PID.m .
clc; clear all; close all
%% Best K parameters
% Best K parameters are as follows:
Kp = 0.005;Kd = 0.001;Ki = 0.007;
Pao_with_PID = Pao_PID(Kp,Kd,Ki);
Pao_no_PID = Pao_PID(); % PID = 1 
N = 1:400;

figure;
plot(N, Pao_with_PID); hold on;
plot(N, Pao_no_PID,'r'); grid on ;
xlim([0 250]);
xlabel('No. of cycles'); ylabel('Average Aortic Pressure [mmHg]');
legend('Avg Aortic Pressure with PID control', 'Avg Aortic Pressure without PID control');
title('Average Aortic Pressure with and without PID control');  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Different K Values %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kp
Kp_p=0.005*3; Kd_p=0.001; Ki_p=0.0001;
G1=Pao_PID(Kp_p,Kd_p,Ki_p);
Kp_p=0.005/2;
G2=Pao_PID(Kp_p,Kd_p,Ki_p);
Kp_p=0.005/20;
G3=Pao_PID(Kp_p,Kd_p,Ki_p);

figure; 
plot(N, Pao_with_PID); hold on
plot(N, G1,'m');plot(N, G2,'g');
plot(N, G3,'r');grid on ;
xlim([0 250]);
xlabel('No. of cycles'); ylabel('Average Aortic Pressure [mmHg]');
legend('Chosen K', 'Kp = 0.015', 'Kp = 0.0025', 'Kp = 0.00025');
title('Average Aortic Pressure: Different Kp Values');  hold off
 
%% Ki
Kp_p=0.005;Kd_p=0.001;Ki_p=0.001;
G1=Pao_PID(Kp_p,Kd_p,Ki_p);
Ki_p=0.00001;
G2=Pao_PID(Kp_p,Kd_p,Ki_p);
Ki_p=0.003;
G3=Pao_PID(Kp_p,Kd_p,Ki_p);

figure;
plot(Pao_with_PID); hold on
plot(G1,'m');plot(G2,'g');
plot(N, G3,'r');grid on ;
xlim([0 250]);
xlabel('No. of cycles'); ylabel('Average Aortic Pressure [mmHg]');
legend('Chosen K', 'Ki = 0.001', 'Ki = 0.00001', 'Ki=0.003');
title('Average Aortic Pressure: Different Ki Values'); hold off

%% Kd
Kp_p=0.005;Kd_p=0.001*7;Ki_p=0.0001;
G1=Pao_PID(Kp_p,Kd_p,Ki_p);
Kd_p=0.001*10;
G2=Pao_PID(Kp_p,Kd_p,Ki_p);
Kd_p=0.0001;
G3=Pao_PID(Kp_p,Kd_p,Ki_p);

figure;
plot(Pao_with_PID); hold on
plot(G1,'m');plot(G2,'g');
plot(N, G3,'r');grid on ;
xlim([0 250]);
xlabel('No. of cycles'); ylabel('Average Aortic Pressure [mmHg]');
legend('Chosen K', 'Kd = 0.007', 'Kd = 0.01','Kd=0.0001');
title('Average Aortic Pressure: Different Kd Values'); hold off
