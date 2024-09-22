function [n_dt,n_simulations]=Martingale_Test(rates,F0,sigmaEq,sigmaPr,T,...
                              seed)
% Function to compute the Martingale Test,
%
% INPUTS
% rates :           Vector of discrete annual compuended rates
% F0 :              Initial Fund value
% sigmeEq :         Equity's volatility
% sigmaPr :         Property's volatility
% T :               Time horizon
% seed :            Seed for random generation
%
% OUTPUTS
% n_dt :            Best number of time steps
% n_simulations :   Best number of simulations

rates=[0;rates];
rates = rates(1:T+1);
% setting a time of simulations and a number of steps
N_sim=1e4:1e4:1e5;
N_timesteps=50:50:200;
err=zeros(length(N_sim),length(N_timesteps));
for i=1:length(N_sim)    
    Number_Simulations=N_sim(i);
    for j=1:length(N_timesteps)
        Number_TimeSteps=N_timesteps(j);
        % scenarios generation
        rng(seed) 
        g=randn(Number_Simulations,Number_TimeSteps);
        time_step=linspace(0,T,Number_TimeSteps+1)';
        length_dt= T/Number_TimeSteps;
        % interpolating the rates on the dates we just found 
        rates_interp=interp1((0:T)',rates,time_step);
        discounts = (1+rates_interp).^(-(1:length(rates_interp))');
        % Discounts
        fwd_discounts=(discounts(2:end)./discounts(1:end-1))';
        fwd_rates=-log(fwd_discounts);
        S_eq=zeros(Number_Simulations,Number_TimeSteps+1);
        S_pr=zeros(Number_Simulations,Number_TimeSteps+1); 
        S_eq(:,1)=F0; % Equity
        S_pr(:,1)=F0; % Property

        for k=2:Number_TimeSteps+1
             S_eq(:,k)=S_eq(:,k-1).*exp((fwd_rates(k-1)-sigmaEq^2/2)*...
                              length_dt+sigmaEq*sqrt(length_dt).*g(:,k-1));
             S_pr(:,k)=S_pr(:,k-1).*exp((fwd_rates(k-1)-sigmaPr^2/2)*...
                              length_dt+sigmaPr*sqrt(length_dt).*g(:,k-1));
        end
        S=0.8*S_eq+0.2*S_pr;
        err(i,j) = norm(-F0*exp(rates_interp.*time_step)+ mean(S)');
        plot(0:length_dt:T,-F0*exp(rates_interp.*time_step)+mean(S)')
        hold on
    end
end

plot(0:T,zeros(T+1,1))
title("Assets simulation")
% hold on 
% legend('n\_simulation = 1e4, n\_time\_step = 1e2',...
%     'n\_simulation = 1e4, n\_time\_step = 1e3',...
%     'n\_simulation = 1e5, n\_time\_step = 1e2',...
%     'n\_simulation = 1e5, n\_time\_step = 1e3','error=0');
% hold off
[minimumValue, linearIndex] = min(err(:));
[rowIndex, colIndex] = ind2sub(size(err), linearIndex);
n_dt=N_timesteps(colIndex);
n_simulations=N_sim(rowIndex);
fprintf(['The minimum error is %d, the best number of time steps is %d and ' ...
    '               of simulations %d.\n'], minimumValue, n_dt, n_simulations);
end 