tic 
clc; clear all; close all;


Data = importdata('V_final_CI.mat');

Nt = 2000; % % number of individuals per trial
deli = 0.2; % % binning log10(NT50) 
Datasave = [0:deli:4]'; % % range of log10(NT50) 

runs = size(Data,2)/Nt; % % number of trials = 5

for i0 = 1:1:runs
c0 = 0;

% % Finding Vpeak and NT50 for each patient
for i = ((i0-1)*Nt)+1:1:i0*Nt
c0 = c0+1;
Pred(c0,:) = [Data{i}.NT50_it max(Data{i}.Vsave(:,2)) Data{i}.Vmax];
end

idp = isnan(Pred(:,1));
Pred(idp,:)=[];

% % Binning and computing protection vs NT50
ct = 0;
for i = 0:deli:4
ct = ct+1;
idx = find(log10(Pred(:,1))>=i & log10(Pred(:,1))<i+deli);
Pred1(ct,:) = [i+(0.5*deli) length(find(Pred(idx,2)>100)) length(idx) (1-(length(find(Pred(idx,2)>100))/length(idx)))*100];
Pred(idx,:)=[];
end

Datasave = [Datasave Pred1(:,4)];

clearvars -except Data Datasave Nt deli i0 runs
end

Datasave(17:end,:)=[];

% % Compute and plot mean and SD
MeanP = mean(Datasave(:,2:runs+1),2)
StdP = std(Datasave(:,2:runs+1)')'

curve1 = MeanP + StdP;
curve2 = MeanP - StdP;

figure
semilogx(10.^Datasave(:,1),MeanP,'-k','linewidth',2)
hold on
plot(10.^Datasave(:,1),curve1,'--k',10.^Datasave(:,1),curve2,'--k','linewidth',1)
xlabel('NT_5_0')
ylabel('Protection (%)')
ylim([-5,105])
set(gca,'linewidth',2,'FontName','Arial','fontsize',12)
axis square

toc
