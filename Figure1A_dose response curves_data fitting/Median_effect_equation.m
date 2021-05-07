function Median_effect_equation

clear all; close all; clc

% % Load data
data = xlsread('Data_BD236.xlsx'); 

% % Remove NaN values (corresponding to fa<1% and fa>99%; see Figure 1 legend)
idx = isnan(data(:,3));
doseresponse=[10.^data(:,1) data(:,3)*100];
doseresponse(idx,:)=[];

% % beta0 provides intial guesses. Note that IC50 is log transformed 
[~,idx]=min(abs(doseresponse(:,2)-50));
init_guess = doseresponse(idx,1);
beta0 = [1 log10(init_guess)]; 

% % parameter and 95% CI estimation 
[betahat,resid,J]=nlinfit(doseresponse(:,1),log10(doseresponse(:,2)./(100-doseresponse(:,2))),@calc,beta0)
betaci = nlparci(betahat,resid,J)

% % Show estimates in command window
m_estimate = betahat(1) 
m_95CI = betaci(1,1:2) 
IC50_estimate = 10^betahat(2) 
IC50_95CI = 10.^betaci(2,1:2)

% % Predict and plot
ln_Dp = [-5:0.1:5]; 
Dp = 10.^ln_Dp; 
Pred = calc(betahat,Dp);

figure
plot(log10(Dp),Pred,'-',log10(doseresponse(:,1)),log10(doseresponse(:,2)./(100-doseresponse(:,2))),'.','MarkerSize',25,'linewidth',2);
xlabel('log_1_0 [BD-236] (\mug ml^-^1)')
ylabel('log_1_0 [f_a/f_u]')
ylim([-3,3])
xlim([-4.3,2.3])
set(gca,'FontSize',18)


function F=calc(beta,D)
m = beta(1);
IC50 = 10^beta(2);
F = m*log10(D) - m*log10(IC50);