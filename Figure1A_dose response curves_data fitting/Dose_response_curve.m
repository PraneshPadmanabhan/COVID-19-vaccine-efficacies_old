function Dose_response_curve

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
[betahat,resid,J]=nlinfit(doseresponse(:,1),doseresponse(:,2),@calc,beta0);
betaci = nlparci(betahat,resid,J);

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
plot(log10(Dp),Pred,'-k',log10(doseresponse(:,1)),doseresponse(:,2),'.r','MarkerSize',35,'linewidth',2);
hold on 
plot(data(1,1),data(1,2)*100,'or',data(2,1),data(2,2)*100,'or',data(9,1),data(9,2)*100,'or','MarkerSize',8,'linewidth',2)
xlabel('log_1_0 [BD-236] (\mug ml^-^1)')
ylabel('Fraction affected, f_a (%)')
ylim([-20,120])
xlim([-4.3,2.3])
set(gca,'FontSize',18)


function Fa=calc(beta,D)
m = beta(1);
IC50 = 10^beta(2);
Fa = (D.^m ./ (IC50^m + D.^m))*100;