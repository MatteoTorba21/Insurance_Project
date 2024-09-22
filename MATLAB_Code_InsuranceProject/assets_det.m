function [Eq,Pr,S_comm,S_profit] = assets_det(Eq0,Pr0,T,rates,RD,COMM)
% Function to simulate the stochastic dynamics of the equity and of the
% property via a GBM.
% 
% INPUTS:
% Eq0:      Equity's value in t0
% Pr0:      Property's value in t0
% T:        Number of years
% rates:    Zero rates curve
% RD:       Regular deduction from the fund value
% COMM:     Commisions to be computed on the fund value    
%
% OUTPUTS:
% Eq:       Matrix of simulated equity values
% Pr:       Matrix of simulated property values
% S_comm:   Commissions matrix: element in position i,j represents the
%           commission in the i-th simulation and in the j-th year
% S_profit: Profits matrix: element in position i,j represents the
%           profit in the i-th simulation and in the j-th year

% Initialize the output variables:
Eq = zeros(1,T+1);
Pr = zeros(1,T+1);
Eq_comm = zeros(1,T+1);
Pr_comm = zeros(1,T+1);
Eq_profit = zeros(1,T+1);
Pr_profit = zeros(1,T+1);
Eq(:,1) = Eq0;
Pr(:,1) = Pr0;

for i=1:T
    % Equity and property values at the i-th year for every simulation:
    Eq(:,i+1) = Eq(:,i)*exp(rates(i))*(1-RD);
    Pr(:,i+1) = Pr(:,i)*exp(rates(i))*(1-RD);
    % Commissions on equity and property values at the i-th year for every simulation:
    Eq_comm(:,i+1) = Eq(:,i)*exp(rates(i))*COMM;
    Pr_comm(:,i+1) = Pr(:,i)*exp(rates(i))*COMM;
    % Compute the profit at the i-th year for every simulation:
    Eq_profit(:,i+1) = Eq(:,i)*exp(rates(i))*(RD-COMM);
    Pr_profit(:,i+1) = Pr(:,i)*exp(rates(i))*(RD-COMM);
end

% Matrix of total commissions computed on the total fund value:
S_comm = Eq_comm(:,2:end) + Pr_comm(:,2:end);

% Matrix of profits:
S_profit = Eq_profit(:,2:end) + Pr_profit(:,2:end);

end