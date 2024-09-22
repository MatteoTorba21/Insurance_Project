function [discounts, fwd_discounts, fwd_rates] = ...
                                        Compute_Df_Fwddf_Fwdrates(rates)
% Function to compute the discounts, the forward discounts and the forward
% rates from the spot rates.
%
% INPUT:
% rates:            Vector of discrete annual compuended rates
%
% OUTPUTS:
% discounts:        Vector of annual discounts
% fwd_discounts:    Vector of annual forward discounts
% fwd_rates:        Vector of forward rates computed in t0

% Compute the discounts:
% discounts = exp(-rates.*(1:length(rates))');
discounts = (1+rates).^(-(1:length(rates))');
% Compute the forward discounts:
fwd_discounts = discounts./[1;discounts(1:end-1)];
% Compute the forward rates:
fwd_rates = -log(fwd_discounts);

end