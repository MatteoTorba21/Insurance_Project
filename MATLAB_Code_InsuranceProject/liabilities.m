function [Liab_death, Liab_lapse, Liab_survive, Expense, Commissions, ...
                  Profits] = liabilities(S,T,lapse,q,x,Pen,C0,infl,expen,...
                  Comm_mat,Profit_mat)

% Function to compute the liabilities and the profits.
%
% INPUTS:
% S:            Assets value
% T:            Time (in years) to the end of the contract
% lapse:        Lapse annual probability vector
% q:            Death probability
% x:            Age of the insured person
% Pen:          Penalties to pay in case of lapse of the insured
% C0:           Invested premium value
% infl:         Yearly inflation rate
% expen:        Yearly expenses
% Comm_mat:     Commissions matrix: element in position i,j represents the
%               commission in the i-th simulation and in the j-th year
% Profit_mat:   Profits matrix: element in position i,j represents the
%               profit in the i-th simulation and in the j-th year
% 
% OUTPUTS:
% Liab_death:   Vector of yearly liabilities caused by death
% Liab_lapse:   Vector of yearly liabilities caused by lapse
% Liab_survive: Vector of yearly liabilities caused by survival and not lapse
%               of the insured until the end of the contract
% Expense:      Vector of yearly liabilities caused by expenses
% Commissions:  Vector of yearly liabilities caused by commissions
% Profits:      Vector of profits

% Initialization of the survival and the "not lapse" probabilities:
surv_prob = zeros(T+1,1); % probability to survive from x to x+t for each t
surv_prob(1) = 1;
prob_not_lapse = zeros(T+1,1); % probability to survive from x to x+t for 
                               % each t
prob_not_lapse(1) = 1;

% Compute the survival probabilities and the "not lapse" probabilities:
for i=1:T
    surv_prob(i+1) = surv_prob(i)*(1-q(x+i));
    prob_not_lapse(i+1) = prob_not_lapse(i)*(1-lapse(i)); % Consider the 
                                                          % cumulative product
                                                          % of the non-lapse
                                                          % probability
end

% Compute lapse and death benefits
[C_lapse,C_death] = funds(S,Pen,C0);

% Liabilities in case of death:
Liab_death = mean(C_death.*(q((x+1):(x+T)).*(prob_not_lapse(1:end-1).*...
                surv_prob(1:end-1)))',1);
% Liabilities in case of lapse:
Liab_lapse = mean(C_lapse.*(lapse'.*(prob_not_lapse(1:end-1).*surv_prob...
                  (2:end)))',1);

% Liabilities in case of survival:
Liab_survive = mean(surv_prob(end)*S(:,end)*prob_not_lapse(end),1);
% Liabilities due to expenses:
Expense = (expen*(1+infl).^(0:T-1)'.*prob_not_lapse(2:end).*surv_prob(2:end))';
% Liabilities due to commissions:
Commissions = mean(Comm_mat.*(prob_not_lapse(2:end).*surv_prob(2:end))',1);

% Profits:
Profits = mean(Profit_mat.*(prob_not_lapse(2:end).*surv_prob(2:end))',1);

end