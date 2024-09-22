%% ANTITHETIC VARIABLES SIMUILATION
% The function BOF is called everytime with the flag variable of the asset
% simlation equal to 2 in order to simulate the assets via the antithetic
% variables technique
clc, clear, close all;
format long
seed = 10;

%% Read data and compute the discounts
% Save life tables' data in a struct:
% MaleData = readLifeTables('Tavole.xlsx',1); % Males life tables
% AllData = readLifeTables('Tavole.xlsx',2);  % Males and females life tables
load("MaleData.mat")
load("AllData.mat")

% Save rates' data:
% rates = readRatesData('EIOPA_RFR_20240331_Term_Structures.xlsx');
load("rates.mat")

%%
% Compute the discounts, the forward discounts and the forward rates:
[discounts, fwd_discounts, fwd_rates] = Compute_Df_Fwddf_Fwdrates...
                                        (rates.spot);


%% Data
% Assets
F0 = 1e+5;     % fund value
C0 = F0;       % invested premium value
Eq0 = 0.8*F0;  % Equity in t0
Pr0 = 0.2*F0;  % Property in t0
sigmaEq = 0.2; % Equity's volatility in GBM dynamics
sigmaPr = 0.1; % Property's volatility in GBM dynamics

x = 60; % Age of the insured
T = 50; % Years

% Liabilities
Pen = 20;         % Penalty in case of lapes
RD = 0.022;       % Regular deduction
COMM = 0.014;     % Commisions percentage
lapse = 0.15*ones(1,T); % Lapse annual probability vector
expenses = 50;    % Yearly expenses
inflation = 0.02; % Yearly inflation

%% Probabilities
q = MaleData.qx/1000; % death probability
p = 1-q;              % survival probability

%% Martingale test:
[opt_t,opt_sim] = Martingale_Test(rates.spot,F0,sigmaEq,sigmaPr,T,seed);
figure
[opt_t_1,opt_sim_1] = Martingale_Test_1(rates.spot,F0,sigmaEq,sigmaPr,...
                        T,seed);
% Set the number of simulations equal to the otpimal value found in the
% Martingale test:
M = opt_sim_1;

%% Base scenario
% Base is a struct which contains all the liabilities and the BOF
rng(seed)
% Simulate assets, compute liabilities and BOF:
base = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation,expenses,F0,discounts,2);

% Esimation of the LEAK:
LEAK = F0-base.liab-base.Profit;
% Proxy of the profit:
profit_proxy = (RD-COMM)*(F0)*base.Duration;
% rates_for_proxy = interp1(1:length(rates.spot),rates.spot,1:floor...
%                    (base.Duration));
% profit_proxy = sum((RD-COMM)*F0*(1+rates_for_proxy).^(-(1:floor(...
%                    base.Duration))));


% Display the results:
fprintf('\nBase scenario data:\n')
fprintf('- Liabilities: %.8f\n',base.liab)
fprintf('- Duration: %.8f\n',base.Duration)
fprintf('- BOF: %.8f\n',base.BOF)
fprintf('- PVFP: %.8f\n',base.Profit)
fprintf('- PVFP proxy: %.8f\n',profit_proxy)
fprintf('- Estimated LEAK: %.8f\n\n',LEAK)

%% Equity stress - Type 1
% Decrese in equity type 1:
s1 = 0.39; %+ 0.0525;
rng(seed)
% Simulate assets, compute liabilities and BOF in the equity stress type 1:
eq1 = BOF(Eq0*(1-s1),Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse,q,x,Pen,Eq0*(1-s1)+Pr0,inflation,expenses,...
                  Eq0*(1-s1)+Pr0,discounts,2);
% Compute the SCR in the equity stress type 1:
SCR_s1 = max(base.BOF-eq1.BOF,0);

% Display the results:
disp('Equity stress type 1 data:')
fprintf('- Liabilities: %.8f\n',eq1.liab)
fprintf('- Duration: %.8f\n',eq1.Duration)
fprintf('- BOF: %.8f\n',eq1.BOF)
fprintf('- SCR: %.8f\n\n',SCR_s1)

%% Equity stress - Type 2
% Decrese in equity type 2:
s2 = 0.49;
rng(seed)
% Simulate assets, compute liabilities and BOF in the equity stress type 2:
eq2 = BOF(Eq0*(1-s2),Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse,q,x,Pen,Eq0*(1-s2)+Pr0,inflation,expenses,...
                  Eq0*(1-s2)+Pr0,discounts,2)
% Compute the SCR in the equity stress type 2:
SCR_s2 = max(base.BOF-eq2.BOF,0);

Mat = [1, 0.75;
       0.75, 1];
% SCR equity:
SCR_eq = sqrt([SCR_s1, SCR_s2]*Mat*[SCR_s1, SCR_s2]');

%% Property stress
% Property shock:
P_sp = 0.25;
rng(seed)
% Simulate assets, compute liabilities and BOF in the property shock case:
pr = BOF(Eq0,Pr0*(1-P_sp),M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse,q,x,Pen,Eq0+Pr0*(1-P_sp),inflation,...
                  expenses,Eq0+Pr0*(1-P_sp),discounts,2)
% Compute the SCR:
SCR_sp = max(base.BOF-pr.BOF,0);

%% Interest rates UP
% Shift of the interest rates taken from EIOPA rates' tables
% Compute the discounts, the forward discounts and the forward rates from 
% the shifted rates:
[discounts_irup, fwd_discounts_irup, fwd_rates_irup] = ...
                                    Compute_Df_Fwddf_Fwdrates(rates.shockup);
rng(seed)
% Simulate assets, compute liabilities and BOF in the IR shock up case:
I_sup = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_irup,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation,expenses,F0,...
                  discounts_irup,2)
% Compute the SCR:
SCR_irup = max(base.BOF-I_sup.BOF,0);

%% Interest rates DOWN
% Shift of the interest rates taken from EIOPA rates' tables
% Compute the discounts, the forward discounts and the forward rates 
% from the shifted rates:
[discounts_irdown, fwd_discounts_down, fwd_rates_down] = ...
                            Compute_Df_Fwddf_Fwdrates(rates.shockdown);
rng(seed)
% Simulate assets, compute liabilities and BOF in the IR shock down case:
I_down = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_down,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation,expenses,...
                  F0,discounts_irdown,2)
% Compute the SCR:
SCR_irdown = max(base.BOF-I_down.BOF,0);

%% EXPENSE stress
rng(seed) % set the seed
% Simulate assets, compute liabilities and BOF after shifting inflation and
% expenses:
exp_up = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation+0.01,expenses*1.1,F0,...
                  discounts,2)
% Compute the SCR:
SCR_exp = max(base.BOF-exp_up.BOF,0);

%% Mortality stress
rng(seed) % set the seed
% Simulate assets, compute liabilities and BOF after shifting the
% mortality:
mort = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse,q*1.15,x,Pen,C0,inflation,expenses,F0,...
                  discounts,2)
% Compute the SCR:
SCR_sm = max(base.BOF-mort.BOF,0);

%% Lapse stress UP
rng(seed) % set the seed
% Lapse stressed up:
lapse_lu = min(lapse*1.5,1);
% Simulate assets, compute liabilities and BOF with the stressed Lapse:
lapse_up = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse_lu,q,x,Pen,C0,inflation,expenses,F0,...
                  discounts,2)
% Compute the SCR:
SCR_lu = max(base.BOF-lapse_up.BOF,0);

%% Lapse stress DOWN
rng(seed) % set the seed
% Lapse stressed down:
lapse_ld = max(lapse*0.5,lapse-0.2)
% Simulate assets, compute liabilities and BOF with the stressed Lapse:
lapse_dw = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse_ld,q,x,Pen,C0,inflation,expenses,F0,...
                  discounts,2)
% Compute the SCR:
SCR_ld = max(base.BOF-lapse_dw.BOF,0);

%% Lapse stress MASS
rng(seed)   % set the seed
% Stress Lapse MASS:
lapse_mass = zeros(size(lapse));
lapse_mass(1) = lapse(1)+0.4;
lapse_mass(2:end) = lapse(2:end);
% Simulate assets, compute liabilities and BOF with the stressed Lapse:
lapse_ms = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse_mass,q,x,Pen,C0,inflation,expenses,F0,...
                  discounts,2)
% Compute the SCR of the Lapse MASS shock:
SCR_mass = max(base.BOF-lapse_ms.BOF,0);

% Compute the SCR for the lapse shocks:
SCR_lapse = max(SCR_lu,max(SCR_mass,SCR_ld));

%% CAT scenario
rng(seed)
q1 = q;
q1(x+1) = q(x+1) + 0.0015;
% Simulate assets, compute liabilities and BOF:
CAT = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse,q1,x,Pen,C0,inflation,expenses,F0,...
                  discounts,2)
% Compute the SCR:
SCR_CAT = max(base.BOF-CAT.BOF,0);

%% RM MKT
A = 0.5; % 0 UP, 0.5 DOWN
CorrM = [1, A, A; % IR
        A, 1, 0.75; % Eq
        A, 0.75, 1]; % Pr

SCR_mkt = sqrt([SCR_irdown, SCR_eq, SCR_sp]*CorrM*[SCR_irdown, SCR_eq,...
                SCR_sp]');

% A = 0; % 0 UP, 0.5 DOWN
% CorrM = [1, A, A; % IR
%         A, 1, 0.75; % Eq
%         A, 0.75, 1]; % Pr
% 
% SCR_mkt = sqrt([SCR_irup, SCR_eq, SCR_sp]*CorrM*[SCR_irup, SCR_eq, SCR_sp]');

%% RM LIFE
CorrL = [1, 0, 0.25, 0.25; % Mort
         0, 1, 0.5, 0.25; % Lap
         0.25, 0.5, 1, 0.25; % Exp
         0.25, 0.25, 0.25, 1]; % CAT

SCR_life = sqrt([SCR_sm, SCR_lapse, SCR_exp, SCR_CAT]*CorrL*[SCR_sm,...
                SCR_lapse, SCR_exp, SCR_CAT]');

%% Matrix Mkt - Life
CorrML = [1, 0.25; % Mkt
          0.25, 1]; % Life

BSCR = sqrt([SCR_mkt, SCR_life]*CorrML*[SCR_mkt, SCR_life]');