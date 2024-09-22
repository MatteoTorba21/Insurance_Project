function [Eq,Pr,S_comm,S_profit] = assets_simulation(Eq0,Pr0,M,T,sigmaEq,sigmaPr,fwd,RD,COMM,flag)
% Function to simulate the stochastic dynamics of the equity and of the
% property via a GBM.
% 
% INPUTS:
% Eq0:      Equity's value in t0
% Pr0:      Property's value in t0
% M:        Number of simulations
% T:        Number of years
% sigmaEq:  Equity's GBM volatility
% sigmaPr:  Property's GBM volatility
% fwd:      Forward rates curve
% RD:       Regular deduction from the fund value
% COMM:     Commisions to be computed on the fund value 
% flag:     Flag variable to choose the technique to do the simulation:
%           -if flag==1, do the 'classic' simulation
%           -if flag==2, do the simulation via antithetic variables technique
%           -if flag==3 do the deterministic projection
%
% OUTPUTS:
% Eq:       Matrix of simulated equity values
% Pr:       Matrix of simulated property values
% S_comm:   Commissions matrix: element in position i,j represents the
%           commission in the i-th simulation and in the j-th year
% S_profit: Profits matrix: element in position i,j represents the
%           profit in the i-th simulation and in the j-th year

% Initialize the output variables:
Eq = zeros(M,T+1);
Pr = zeros(M,T+1);
Eq_comm = zeros(M,T+1);
Pr_comm = zeros(M,T+1);
Eq_profit = zeros(M,T+1);
Pr_profit = zeros(M,T+1);
Eq(:,1) = Eq0;
Pr(:,1) = Pr0;
% Vector of times:
dt = ones(1,T);

switch (flag)
    case 1  % 'Classic' simulation
        % Normal random variables:
        g = randn(M,T);
        g1 = randn(M,T);
        for i=1:T
            % Equity and property values at the i-th year for every simulation:
            Eq(:,i+1) = Eq(:,i).*exp((fwd(i)-sigmaEq^2/2)*dt(i)+sigmaEq*sqrt(dt(i))*g(:,i))*(1-RD);
            Pr(:,i+1) = Pr(:,i).*exp((fwd(i)-sigmaPr^2/2)*dt(i)+sigmaPr*sqrt(dt(i))*g1(:,i))*(1-RD);
            % Commissions on equity and property values at the i-th year for every simulation:
            Eq_comm(:,i+1) = Eq(:,i).*exp((fwd(i)-sigmaEq^2/2)*dt(i)+sigmaEq*sqrt(dt(i))*g(:,i))*COMM;
            Pr_comm(:,i+1) = Pr(:,i).*exp((fwd(i)-sigmaPr^2/2)*dt(i)+sigmaPr*sqrt(dt(i))*g1(:,i))*COMM;
            % Compute the profit at the i-th year for every simulation:
            Eq_profit(:,i+1) = Eq(:,i).*exp((fwd(i)-sigmaEq^2/2)*dt(i)+sigmaEq*sqrt(dt(i))*g(:,i))*(RD-COMM);
            Pr_profit(:,i+1) = Pr(:,i).*exp((fwd(i)-sigmaPr^2/2)*dt(i)+sigmaPr*sqrt(dt(i))*g1(:,i))*(RD-COMM);

        end
    
    case 2  % Simulation via antithetic variables technique
        % Normal random variables:
        g = randn(M/2,T);
        g1 = randn(M/2,T);

        for i=1:T
            % Equity and property values at the i-th year for every simulation:
            Eq(1:M/2,i+1) = Eq(:,i).*exp((fwd(i)-sigmaEq^2/2)*dt(i)+sigmaEq*sqrt(dt(i))*g(:,i))*(1-RD);
            Pr(1:M/2,i+1) = Pr(:,i).*exp((fwd(i)-sigmaPr^2/2)*dt(i)+sigmaPr*sqrt(dt(i))*g1(:,i))*(1-RD);
            Eq(M/2+1:M,i+1) = Eq(:,i).*exp((fwd(i)-sigmaEq^2/2)*dt(i)-sigmaEq*sqrt(dt(i))*g(:,i))*(1-RD);
            Pr(M/2+1:M,i+1) = Pr(:,i).*exp((fwd(i)-sigmaPr^2/2)*dt(i)-sigmaPr*sqrt(dt(i))*g1(:,i))*(1-RD);
            % Commissions on equity and property values at the i-th year for every simulation:
            Eq_comm(1:M/2,i+1) = Eq(:,i).*exp((fwd(i)-sigmaEq^2/2)*dt(i)+sigmaEq*sqrt(dt(i))*g(:,i))*COMM;
            Pr_comm(1:M/2,i+1) = Pr(:,i).*exp((fwd(i)-sigmaPr^2/2)*dt(i)+sigmaPr*sqrt(dt(i))*g1(:,i))*COMM;
            Eq_comm(M/2+1:M,i+1) = Eq(:,i).*exp((fwd(i)-sigmaEq^2/2)*dt(i)-sigmaEq*sqrt(dt(i))*g(:,i))*COMM;
            Pr_comm(M/2+1:M,i+1) = Pr(:,i).*exp((fwd(i)-sigmaPr^2/2)*dt(i)-sigmaPr*sqrt(dt(i))*g1(:,i))*COMM;
            % Compute the profit at the i-th year for every simulation:
            Eq_profit(1:M/2,i+1) = Eq(:,i).*exp((fwd(i)-sigmaEq^2/2)*dt(i)+sigmaEq*sqrt(dt(i))*g(:,i))*(RD-COMM);
            Pr_profit(1:M/2,i+1) = Pr(:,i).*exp((fwd(i)-sigmaPr^2/2)*dt(i)+sigmaPr*sqrt(dt(i))*g1(:,i))*(RD-COMM);
            Eq_profit(M/2+1:M,i+1) = Eq(:,i).*exp((fwd(i)-sigmaEq^2/2)*dt(i)-sigmaEq*sqrt(dt(i))*g(:,i))*(RD-COMM);
            Pr_profit(M/2+1:M,i+1) = Pr(:,i).*exp((fwd(i)-sigmaPr^2/2)*dt(i)-sigmaPr*sqrt(dt(i))*g1(:,i))*(RD-COMM);
        end

    case 3  % Deterministic projection
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
    otherwise
end


% Matrix of total commissions computed on the total fund value:
S_comm = Eq_comm(:,2:end) + Pr_comm(:,2:end);

% Matrix of profits:
S_profit = Eq_profit(:,2:end) + Pr_profit(:,2:end);

end