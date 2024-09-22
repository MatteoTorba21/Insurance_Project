%% Insurance project AY2023-2024
% Callini Matteo, Carzaniga Carlotta, Gaspari Cecilia, Torba Matteo
clear, clc, close all;
format long
seed = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QUESTION 1 AND 2: Compute the Basic Solvency Capital Requirement 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read data and compute the discounts
% Save life tables' data in a struct:
% MaleData = readLifeTables('Tavole.xlsx',1); % Males life tables
% AllData = readLifeTables('Tavole.xlsx',2);  % Males and females life tables
load("MaleData.mat")
load("AllData.mat")

% Save rates' data:
% rates = readRatesData('EIOPA_RFR_20240331_Term_Structures.xlsx');
load("rates.mat")

% Compute the discounts, the forward discounts and the forward rates:
[discounts, fwd_discounts, fwd_rates] = Compute_Df_Fwddf_Fwdrates(rates.spot);

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
inflation = 0.02; % Yearly inflation rate

%% Probabilities
q = MaleData.qx/1000; % death probability
p = 1-q;              % survival probability

%% Martingale Test
tic
[opt_t,opt_sim] = Martingale_Test(rates.spot,F0,sigmaEq,sigmaPr,T,seed);
figure
[opt_t_1,opt_sim_1] = Martingale_Test_1(rates.spot,F0,sigmaEq,sigmaPr,T,seed);
% Set the number of simulations equal to the otpimal value found in the
% Martingale test:
M = opt_sim_1;
fprintf("Martingale test computation time: ")
toc

%% Base scenario
% Base is a struct which contains all the liabilities and the BOF
rng(seed)
% Simulate assets, compute liabilities and BOF:
base = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation,expenses,F0,discounts);
% Display the results:
fprintf('\n --- Point a and Point b ---\n\n')
disp('Base scenario data:')
fprintf('- Liabilities: %.8f\n',base.liab)
fprintf('- Duration: %.8f\n',base.Duration)
fprintf('- BOF: %.8f\n\n',base.BOF)

%% Equity stress - Type 1
% Decrese in equity type 1
s1 = 0.3950 + 0.0525; % We consider the symmetric adjustment term (5.25%)
rng(seed)
% Simulate assets, compute liabilities and BOF in the equity stress type 1:
eq1 = BOF(Eq0*(1-s1),Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse,q,x,Pen,Eq0*(1-s1)+Pr0,inflation,expenses,...
                  Eq0*(1-s1)+Pr0,discounts);
% Compute the SCR in the equity stress type 1:
SCR_s1 = max(base.BOF-eq1.BOF,0);

% Display the results:
disp('Equity stress type 1 data:')
fprintf('- Liabilities: %.8f\n',eq1.liab)
fprintf('- Duration: %.8f\n',eq1.Duration)
fprintf('- BOF: %.8f\n',eq1.BOF)
fprintf('- SCR: %.8f\n\n',SCR_s1)

%% Property stress
% Property shock:
P_sp = 0.25;
rng(seed)
% Simulate assets, compute liabilities and BOF in the property shock case:
pr = BOF(Eq0,Pr0*(1-P_sp),M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse,q,x,Pen,Eq0+Pr0*(1-P_sp),inflation,expenses,...
                  Eq0+Pr0*(1-P_sp),discounts);
% Compute the SCR:
SCR_sp = max(base.BOF-pr.BOF,0);

% Display the results:
disp('Property stress data:')
fprintf('- Liabilities: %.8f\n',pr.liab)
fprintf('- Duration: %.8f\n',pr.Duration)
fprintf('- BOF: %.8f\n',pr.BOF)
fprintf('- SCR: %.8f\n\n',SCR_sp)

%% Interest rates UP
% Shift of the interest rates taken from EIOPA rates' tables
% Compute the discounts, the forward discounts and the forward rates from 
% the shifted rates:
[discounts_irup, fwd_discounts_irup, fwd_rates_irup] = Compute_Df_Fwddf_Fwdrates(...
                                                        rates.shockup);
rng(seed)
% Simulate assets, compute liabilities and BOF in the IR shock up case:
I_sup = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_irup,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation,expenses,F0,discounts_irup);
% Compute the SCR:
SCR_irup = max(base.BOF-I_sup.BOF,0);

% Display the results:
disp('Interest rates up stress data:')
fprintf('- Liabilities: %.8f\n',I_sup.liab)
fprintf('- Duration: %.8f\n',I_sup.Duration)
fprintf('- BOF: %.8f\n',I_sup.BOF)
fprintf('- SCR: %.8f\n\n',SCR_irup)

%% Interest rates DOWN
% Shift of the interest rates taken from EIOPA rates' tables
% Compute the discounts, the forward discounts and the forward rates from 
% the shifted rates:
[discounts_irdown, fwd_discounts_down, fwd_rates_down] = Compute_Df_Fwddf_Fwdrates(...
                                                         rates.shockdown);
rng(seed)
% Simulate assets, compute liabilities and BOF in the IR shock down case:
I_down = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_down,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation,expenses,F0,discounts_irdown);
% Compute the SCR:
SCR_irdown = max(base.BOF-I_down.BOF,0);

% Display the results:
disp('Interest rates down stress data:')
fprintf('- Liabilities: %.8f\n',I_down.liab)
fprintf('- Duration: %.8f\n',I_down.Duration)
fprintf('- BOF: %.8f\n',I_down.BOF)
fprintf('- SCR: %.8f\n\n',SCR_irdown)

%% EXPENSE stress
rng(seed) % set the seed
% Simulate assets, compute liabilities and BOF after shifting inflation and
% expenses:
exp_up = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation+0.01,expenses*1.1,F0,...
                  discounts);
% Compute the SCR:
SCR_exp = max(base.BOF-exp_up.BOF,0);

% Display the results:
disp('Expense stress data:')
fprintf('- Liabilities: %.8f\n',exp_up.liab)
fprintf('- Duration: %.8f\n',exp_up.Duration)
fprintf('- BOF: %.8f\n',exp_up.BOF)
fprintf('- SCR: %.8f\n\n',SCR_exp)

%% Mortality stress
rng(seed) % set the seed
% Simulate assets, compute liabilities and BOF after shifting the
% mortality:
mort = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse,q*1.15,x,Pen,C0,inflation,expenses,F0,...
                  discounts);
% Compute the SCR:
SCR_sm = max(base.BOF-mort.BOF,0);

% Display the results:
disp('Mortality stress data:')
fprintf('- Liabilities: %.8f\n',mort.liab)
fprintf('- Duration: %.8f\n',mort.Duration)
fprintf('- BOF: %.8f\n',mort.BOF)
fprintf('- SCR: %.8f\n\n',SCR_sm)

%% Lapse stress UP
rng(seed) % set the seed
% Lapse stressed up:
lapse_lu = min(lapse*1.5,1);
% Simulate assets, compute liabilities and BOF with the stressed Lapse:
lapse_up = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse_lu,q,x,Pen,C0,inflation,expenses,F0,discounts);
% Compute the SCR:
SCR_lu = max(base.BOF-lapse_up.BOF,0);

% Display the results:
disp('Lapse stress up data:')
fprintf('- Liabilities: %.8f\n',lapse_up.liab)
fprintf('- Duration: %.8f\n',lapse_up.Duration)
fprintf('- BOF: %.8f\n',lapse_up.BOF)
fprintf('- SCR: %.8f\n\n',SCR_lu)

%% Lapse stress DOWN
rng(seed) % set the seed
% Lapse stressed down:
lapse_ld = max(lapse*0.5,lapse-0.2);
% Simulate assets, compute liabilities and BOF with the stressed Lapse:
lapse_dw = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse_ld,q,x,Pen,C0,inflation,expenses,F0,discounts);
% Compute the SCR:
SCR_ld = max(base.BOF-lapse_dw.BOF,0);

% Display the results:
disp('Lapse stress down data:')
fprintf('- Liabilities: %.8f\n',lapse_dw.liab)
fprintf('- Duration: %.8f\n',lapse_dw.Duration)
fprintf('- BOF: %.8f\n',lapse_dw.BOF)
fprintf('- SCR: %.8f\n\n',SCR_ld)

%% Lapse stress MASS
rng(seed)   % set the seed
% Stress Lapse MASS:
lapse_mass = zeros(size(lapse));
lapse_mass(1) = lapse(1)+0.4;
lapse_mass(2:end) = lapse(2:end);
% Simulate assets, compute liabilities and BOF with the stressed Lapse:
lapse_ms = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse_mass,q,x,Pen,C0,inflation,expenses,F0,...
                  discounts);
% Compute the SCR of the Lapse MASS shock:
SCR_mass = max(base.BOF-lapse_ms.BOF,0);

% Compute the SCR for the lapse shocks:
SCR_lapse = max(SCR_lu,max(SCR_mass,SCR_ld));

% Display the results:
disp('Lapse mass stress data:')
fprintf('- Liabilities: %.8f\n',lapse_ms.liab)
fprintf('- Duration: %.8f\n',lapse_ms.Duration)
fprintf('- BOF: %.8f\n',lapse_ms.BOF)
fprintf('- SCR: %.8f\n\n',SCR_mass)

fprintf('- SCR lapse: %.8f\n\n',SCR_lapse)

%% CAT scenario
rng(seed)
q1 = q;
q1(x+1) = q(x+1) + 0.0015;
% Simulate assets, compute liabilities and BOF:
CAT = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse,q1,x,Pen,C0,inflation,expenses,F0,discounts);
% Compute the SCR:
SCR_CAT = max(base.BOF-CAT.BOF,0);

% Display the results:
disp('CAT stress data:')
fprintf('- Liabilities: %.8f\n',CAT.liab)
fprintf('- Duration: %.8f\n',CAT.Duration)
fprintf('- BOF: %.8f\n',CAT.BOF)
fprintf('- SCR: %.8f\n\n',SCR_CAT)

%% RM MKT
A = 0.5; % 0 UP, 0.5 DOWN
CorrM = [1, A, A; % IR
        A, 1, 0.75; % Eq
        A, 0.75, 1]; % Pr

SCR_mkt = sqrt([SCR_irdown, SCR_s1, SCR_sp]*CorrM*[SCR_irdown, SCR_s1, SCR_sp]');

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

SCR_life = sqrt([SCR_sm, SCR_lapse, SCR_exp, SCR_CAT]*CorrL*[SCR_sm, ...
                SCR_lapse, SCR_exp, SCR_CAT]');

%% Matrix Mkt - Life
CorrML = [1, 0.25; % Mkt
          0.25, 1]; % Life

BSCR = sqrt([SCR_mkt, SCR_life]*CorrML*[SCR_mkt, SCR_life]');

fprintf('- SCR life: %.8f\n\n',SCR_life)
fprintf('- SCR market: %.8f\n\n',SCR_mkt)
fprintf('- BSCR: %.8f\n\n',BSCR)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Point d

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open questions
fprintf('\n --- Point d ---\n\n')
% Save rates' data:
rates_SHIFTED = readRatesData('1bpsSHIFT_UP_EIOPA_RFR_20240331.xlsx');

% Positive shift 
[discounts_up, fwd_discounts_up, fwd_rates_up] = Compute_Df_Fwddf_Fwdrates(...
                                                 rates_SHIFTED.spot);

%% Base scenario
% Set the seed:
rng(seed)
base_up = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_up,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation,expenses,F0,discounts_up,3);

% Display the results:
disp('Base scenario data - Shift +1 bps rates curve')
fprintf('- Liabilities: %.8f\n',base_up.liab)
fprintf('- Duration: %.8f\n',base_up.Duration)
fprintf('- BOF: %.8f\n\n',base_up.BOF)

%% Equity stress - Type 1
rng(seed)
eq1_up = BOF(Eq0*(1-s1),Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_up,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation,expenses,Eq0*(1-s1)+Pr0,...
                  discounts_up);

%% Property stress
rng(seed)
pr_up = BOF(Eq0,Pr0*(1-P_sp),M,T,sigmaEq,sigmaPr,fwd_rates_up,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation,expenses,Eq0+Pr0*(1-P_sp),...
                  discounts_up);

%% Interest rates UP
% Shift of the interest rates taken from EIOPA rates' tables
% Compute the discounts, the forward discounts and the forward rates:
[discounts_irup_su, ~, fwd_rates_irup_su] = Compute_Df_Fwddf_Fwdrates(...
    rates_SHIFTED.shockup);
rng(seed)
I_sup_up = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_irup_su,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation,expenses,F0,...
                  discounts_irup_su);

%% Interest rates DOWN
% Shift of the interest rates taken from EIOPA rates' tables
% Compute the discounts, the forward discounts and the forward rates:
[discounts_irdown_sd, ~, fwd_rates_down_sd] = Compute_Df_Fwddf_Fwdrates(...
    rates_SHIFTED.shockdown);
rng(seed)
I_down_up = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_down_sd,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation,expenses,F0,...
                  discounts_irdown_sd);

%% EXPENSE stress
rng(seed) % set the seed
exp_up_up = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_up,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation+0.01,expenses*1.1,...
                  F0,discounts_up);

%% Mortality stress
rng(seed) % set the seed
mort_up = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_up,...
                  RD,COMM,lapse,q*1.15,x,Pen,C0,inflation,expenses,F0,...
                  discounts_up);

%% Lapse stress UP
rng(seed) % set the seed
lapse_up_up = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_up,...
                  RD,COMM,lapse_lu,q,x,Pen,C0,inflation,expenses,F0,...
                  discounts_up);

%% Lapse stress DOWN
rng(seed) % set the seed
lapse_dw_up = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_up,...
                  RD,COMM,lapse_ld,q,x,Pen,C0,inflation,expenses,F0,...
                  discounts_up);

%% Lapse stress MASS
rng(seed)   % set the seed
lapse_ms_up = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_up,...
                  RD,COMM,lapse_mass,q,x,Pen,C0,inflation,expenses,F0,...
                  discounts_up);

%% CAT scenario
rng(seed)
CAT_up = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_up,...
                  RD,COMM,lapse,q1,x,Pen,C0,inflation,expenses,F0,discounts_up);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Negative shift 
% Save rates' data:
rates_SHIFTED_DOWN = readRatesData('1bpsSHIFT_DOWN_EIOPA_RFR_20240331.xlsx');
[discounts_dw, fwd_discounts_dw, fwd_rates_dw] = Compute_Df_Fwddf_Fwdrates(...
                                                 rates_SHIFTED_DOWN.spot);

%% Base scenario
% Set the seed:
rng(seed)
base_dw = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_dw,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation,expenses,F0,discounts_dw);

% Display the results:
disp('Base scenario data - Shift -1 bps rates curve')
fprintf('- Liabilities: %.8f\n',base_dw.liab)
fprintf('- Duration: %.8f\n',base_dw.Duration)
fprintf('- BOF: %.8f\n\n',base_dw.BOF)

%% Equity stress - Type 1
rng(seed)
eq1_dw = BOF(Eq0*(1-s1),Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_dw,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation,expenses,Eq0*(1-s1)+Pr0,...
                  discounts_dw);

%% Property stress
rng(seed)
pr_dw = BOF(Eq0,Pr0*(1-P_sp),M,T,sigmaEq,sigmaPr,fwd_rates_dw,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation,expenses,Eq0+Pr0*(1-P_sp),...
                  discounts_dw);

%% Interest rates UP
% Shift of the interest rates taken from EIOPA rates' tables
% Compute the discounts, the forward discounts and the forward rates:
[discounts_irup_su, fwd_discounts_irup_su, fwd_rates_irup_su] = ...
                      Compute_Df_Fwddf_Fwdrates(rates_SHIFTED_DOWN.shockup);
rng(seed)
I_sup_dw = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_irup_su,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation,expenses,F0,...
                  discounts_irup_su);

%% Interest rates DOWN
% Shift of the interest rates taken from EIOPA rates' tables
% Compute the discounts, the forward discounts and the forward rates:
[discounts_irdown_sd, fwd_discounts_down_sd, fwd_rates_down_sd] = ...
                    Compute_Df_Fwddf_Fwdrates(rates_SHIFTED_DOWN.shockdown);
rng(seed)
I_down_dw = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_down_sd,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation,expenses,F0,...
                  discounts_irdown_sd);

%% EXPENSE stress
rng(seed) % set the seed
exp_up_dw = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_dw,...
                  RD,COMM,lapse,q,x,Pen,C0,inflation+0.01,expenses*1.1,F0,...
                  discounts_dw);

%% Mortality stress
rng(seed) % set the seed
mort_dw = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_dw,...
                  RD,COMM,lapse,q*1.15,x,Pen,C0,inflation,expenses,F0,...
                  discounts_dw);

%% Lapse stress UP
rng(seed) % set the seed
lapse_up_dw = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_dw,...
                  RD,COMM,lapse_lu,q,x,Pen,C0,inflation,expenses,F0,...
                  discounts_dw);

%% Lapse stress DOWN
rng(seed) % set the seed
lapse_dw_dw = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_dw,...
                  RD,COMM,lapse_ld,q,x,Pen,C0,inflation,expenses,F0,...
                  discounts_dw);

%% Lapse stress MASS
rng(seed)   % set the seed
lapse_ms_dw = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_dw,...
                  RD,COMM,lapse_mass,q,x,Pen,C0,inflation,expenses,F0,...
                  discounts_dw);

%% CAT scenario
rng(seed)
CAT_dw = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_dw,...
                  RD,COMM,lapse,q1,x,Pen,C0,inflation,expenses,F0,...
                  discounts_dw);
%% Plots
% Assets computation
[Eq,Pr,~] = assets(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,RD,COMM);
% Sum the simulated assets:
S = Eq + Pr;

% Assets computation
[Eq_up,Pr_up,~] = assets(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_up,RD,COMM);
% Sum the simulated assets:
S_up = Eq_up + Pr_up;

% Assets computation
[Eq_dw,Pr_dw,~] = assets(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates_dw,RD,COMM);
% Sum the simulated assets:
S_dw = Eq_dw + Pr_dw;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Increasing age
% Base scenario

% Set the seed:
BOF_age = zeros(1,5);
for i=1:5
    rng(seed)
    base_1 = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,...
                      RD,COMM,lapse,q,x+i,Pen,C0,inflation,expenses,...
                      F0,discounts);
    % Display the results:
    fprintf('Base scenario data - Insured age = %d:\n',x+i)
    fprintf('- Liabilities: %.8f\n',base_1.liab)
    fprintf('- Duration: %.8f\n',base_1.Duration)
    fprintf('- BOF: %.8f\n\n',base_1.BOF)
    BOF_age(i) = base_1.BOF;
end

figure
plot(x+1:x+5,BOF_age,'-ro')
title('BOF base case - Increasing age')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Males and females
% 1 man - 1 female 
% Probabilities
q_all = AllData.qx/1000; % death probability
p_all = 1-q_all; % survival probability

figure
plot(0:length(q)-1,q,0:length(q_all)-1,q_all)
legend('Men probabilities','Men and Female probabilities')

%% Base scenario
% Set the seed:
rng(seed)
base_all = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,...
                  RD,COMM,lapse,q_all,x,Pen,C0,inflation,expenses,F0,...
                  discounts);

% Display the results:
disp('Base scenario data - 1 male and 1 female')
fprintf('- Liabilities: %.8f\n',base_all.liab)
fprintf('- Duration: %.8f\n',base_all.Duration)
fprintf('- BOF: %.8f\n\n',base_all.BOF)
