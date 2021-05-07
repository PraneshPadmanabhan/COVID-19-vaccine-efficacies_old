function Viral_dynamics_Fig4B
clear; clc; 
close all

% Efficacies explored------------------------------------------------------ 
epi_CI = [0 0.4 0.6 0.8 0.95 0.99];

% Parameters--------------------------------------------------------------- 
p.B = 10^-(7.22185); 
p.pX = 4;
p.del = 0.6;
p.p = 390; 
p.c = 20;
p.sX = 1;
p.fiX = 100;
p.dX = 0.2;

% Initial conditions-------------------------------------------------------
p.T0 = 3*10^7;
p.I0 = 1;
p.V0 = p.p/p.c;
p.X0 = 0;

% -------------------------------------------------------------------------
tspan = [0:0.01:40];
Datasave = [tspan'];

% % prediction viral dynamics for different efficacies
for i = 1:1:length(epi_CI)
p.epi = epi_CI(i);
N0 = [p.T0 p.I0 p.V0 p.X0];

[T,N] =  ode23s(@(t,y) SARSCoV2(t,y,p),tspan,N0);

Datasave = [Datasave N(:,3)];
assignin('base','V',Datasave);
assignin('base','epi_CI',epi_CI);
% assignin('base','name',N);

hold on
plot(T,log10(N(:,3)))
xlim([0,35]); ylim([0,10])
xlabel('Time (d)'); ylabel('log_1_0 [V]')
set(gca,'fontsize', 18)
end

function dy = SARSCoV2(~,y,p)

dy = zeros(4,1);

T = y(1);
I = y(2);
V = y(3);
X = y(4);

dy(1) = -p.B*(1-p.epi)*T*V - p.pX*X*T;
dy(2) = p.B*(1-p.epi)*T*V - p.del*I;
dy(3) = p.p*I - p.c*V;
dy(4) = p.sX*(I/(p.fiX+I))*(1-X) - p.dX*X;




