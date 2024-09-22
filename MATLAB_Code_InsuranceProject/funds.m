function [C_lapse, C_death] = funds(S,Pen,C0)
% Function to compute the lapse and death benefits of the contract.
%
% INPUTS:
% S:        Assets value
% Pen:      Penalty in case of lapse
% C0:       Initial capital
%
% OUTPUTS:
% C_lapse:  Lapse benefits
% C_death:  Death benefits

% Compute the benefits:
C_lapse = (S(:,2:end)-Pen);
C_death = max(S(:,2:end),C0);

end