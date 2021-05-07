function Model_final_random_IgG_VK_spectrum_CI
clear; clc; close all; tic

% Load data containing NAbs with distinct m and IC50 values 
Data = importdata('m_IC50_ellipse_landscape_v2.mat'); 
IC50_save = 10.^Data(:,1); m_save = Data(:,2); 

% Time post infection for each case
tspan = [0:0.01:30];

Totalrun = 10000;

for i = 1:1:Totalrun

% Parameters with patient-to-patient variations
p.B = 10^-(7.22185+(0.8*(rand-0.5))); 
p.pX = 4+(2*(rand-0.5));
p.del = 0.6+(1*(rand-0.5));
p.p = 390+(200*(rand-0.5)); 
p.c = 20+(10*(rand-0.5));
p.sX = 1+(1*(rand-0.5));
p.fiX = 100+(200*(rand-0.5));
p.dX = 0.2+(0.1*(rand-0.5));    
p.T0 = 3*10^7;
p.I0 = 1;
p.V0 = ((p.p/p.c)*p.I0);
p.X0 = 0;
p.D0 = 10^(0+(2*(rand-0.5)));
p.N = 10;

Gct = i;
IC50 = IC50_save(((Gct-1)*10 + 1):10*Gct,1)'; 
m = m_save(((Gct-1)*10 + 1):10*Gct,1)'; 
Dilution = 1; 

% % Compute concentration of each NAb
Di(1,1:1:p.N) = p.D0/p.N; 

% % Loewe additivity expression
FF=@(epi) [1 - (Di(1)/(Dilution*IC50(1)))*(((1/epi)-1)^(1/m(1))) ...
             - (Di(2)/(Dilution*IC50(2)))*(((1/epi)-1)^(1/m(2)))...
             - (Di(3)/(Dilution*IC50(3)))*(((1/epi)-1)^(1/m(3)))...
             - (Di(4)/(Dilution*IC50(4)))*(((1/epi)-1)^(1/m(4)))...
             - (Di(5)/(Dilution*IC50(5)))*(((1/epi)-1)^(1/m(5)))...
             - (Di(6)/(Dilution*IC50(6)))*(((1/epi)-1)^(1/m(6)))...
             - (Di(7)/(Dilution*IC50(7)))*(((1/epi)-1)^(1/m(7)))...
             - (Di(8)/(Dilution*IC50(8)))*(((1/epi)-1)^(1/m(8)))...
             - (Di(9)/(Dilution*IC50(9)))*(((1/epi)-1)^(1/m(9)))...
             - (Di(10)/(Dilution*IC50(10)))*(((1/epi)-1)^(1/m(10)))];

options = optimoptions('lsqnonlin','Display','none');
[x] = lsqnonlin(FF,[0.1],[0],[1],options);

p.epi = x;
p.IC50_it = IC50;
p.m_it = m;

% Storing NT50 values
if x>0.5
    p.NT50_it = sum(1./IC50)*(p.D0/p.N);
else
    p.NT50_it = NaN;
end

% initial conditions
N0 = [p.T0 p.I0 p.V0 p.X0];

[T,N] =  ode45(@(t,y) SARSCoV2(t,y,p),tspan,N0);

p.Vsave = [tspan' N(:,3)];
p.Vmax = max(N(:,3));

ctt = i;
Datasave{ctt} = p;

clearvars -except Datasave IC50_save m_save tspan Totalrun i

end

save(strcat('V_final_CI.mat'),'Datasave')

toc
end

function dy = SARSCoV2(t,y,p)

dy = zeros(4,1);

T = y(1);
I = y(2);
V = y(3);
M = y(4);

dy(1) = -p.B*(1-p.epi)*T*V - p.pX*M*T;
dy(2) = p.B*(1-p.epi)*T*V - p.del*I;
dy(3) = p.p*I - p.c*V;
dy(4) = p.sX*(I/(p.fiX+I))*(1-M) - p.dX*M;
end
