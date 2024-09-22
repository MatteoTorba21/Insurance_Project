function base = BOF(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd_rates,RD,COMM,lapse,q,...
                    x,Pen,C0,inflation,expenses,F0,discounts,flag)

% Function to simulate the assets, compue the liabilities and the BOF given
% the input parameters.
%
% INPUTS:
% Eq0:          Equity's value in t0
% Pr0:          Property's value in t0
% M:            Number of simulations
% T:            Number of years
% sigmaEq:      Equity's GBM volatility
% sigmaPr:      Property's GBM volatility
% fwd_rates:    Forward rates curve
% RD:           Regular deduction from the fund value
% COMM:         Commisions to be computed on the fund value
% lapse:        Annual lapse probability vector
% q:            Death probability
% x:            Age of the insured person
% Pen:          Penalties applied in case of lapse
% C0:           Invested premium value
% inflation:    Inflation annually rate
% expenses:     Yearly expenses
% F0:           Fund initial value
% discounts:    Vector of discount factors
% flag:         Flag variable to choose the technique to do the simulation:
%               -if flag==1, do the 'classic' simulation
%               -if flag==2, do the simulation via antithetic variables
%                            technique
%               -if flag==3 do the deterministic projection
%
% OUTPUT:
% base:         Struct of the output with the following fields:
%               -base.Liab_death -> vector of yearly liabilities in case of
%                                   death (not discounted)
%               -base.Liab_lapse -> vector of yearly liabilities in case of
%                                   lapse (not discounted)
%               -base.Liab_survive -> vector of yearly liabilities in case of
%                                     no death nor lapse (not discounted)
%               -base.Expense -> vector of yearly expenses (not discounted)
%               -base.Commissions -> vector of yearly commissions 
%                                    (not discounted)

% Set the flag variable such that the 'classic' simulation is done if it is
% not specified differently:
if nargin<19
    flag = 1;
end

% Simulate equity, property and commisions matrix:
if flag == 2 % simulation via antithetic variables technique
    [Eq,Pr,Comm_mat,Profit_mat] = assets_antithetic_variables(Eq0,Pr0,M,T,...
                                    sigmaEq,sigmaPr,fwd_rates,RD,COMM);
elseif flag==3 % deterministic projection
    [Eq,Pr,Comm_mat,Profit_mat] = assets_det(Eq0,Pr0,T,fwd_rates,RD,COMM);
else    % 'classic' simulation
    [Eq,Pr,Comm_mat,Profit_mat] = assets(Eq0,Pr0,M,T,sigmaEq,sigmaPr,...
                                    fwd_rates,RD,COMM);
end
% Sum the simulated assets:
S = Eq + Pr;

% Compute the liabilities:
[base.Liab_death, base.Liab_lapse, base.Liab_survive, base.Expense,...
        base.Commissions,base.Profit] = liabilities(S,T,lapse,q,x,...
                             Pen,C0,inflation,expenses,Comm_mat,Profit_mat);
% Compute the duration:
liab_duration = (((base.Liab_death + base.Liab_lapse + base.Expense +...
                 base.Commissions) + [zeros(1,T-1),base.Liab_survive]).*...
                 (1:T))*discounts(1:T);
disc_liab = ((base.Liab_death + base.Liab_lapse + ...
                base.Expense + base.Commissions) +...
                [zeros(1,T-1),base.Liab_survive])*discounts(1:T);
base.Duration = liab_duration/disc_liab;

% Discount and sum up all the liabilities:
base.Liab_death = base.Liab_death*discounts(1:T);
base.Liab_lapse = base.Liab_lapse*discounts(1:T);
base.Expense = base.Expense*discounts(1:T);
base.Commissions = base.Commissions*discounts(1:T);
base.Liab_survive = base.Liab_survive*discounts(T);
base.liab = (base.Liab_death + base.Liab_lapse + base.Expense + ...
             base.Commissions) + base.Liab_survive;

% Compute the BOF:
base.BOF = F0-base.liab;

% Discount the profits:
base.Profit = base.Profit*discounts(1:T);

end