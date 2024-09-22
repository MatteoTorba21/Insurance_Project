function rates = readRatesData(filename)
% Function to extract data about the rates from EIOPA saved in a .xls file
%
% INPUTS: 
% filename: Name of the file
%
% OUTPUT:
% rates:    Struct with extracted data

% Extract European spot rates (no VA):
rates.spot = xlsread(filename,3,'C11:C160');
% Extract European spot rates (no VA) shocked up:
rates.shockup = xlsread(filename,5,'C11:C160');
% Extract European spot rates (no VA) shocked down:
rates.shockdown = xlsread(filename,6,'C11:C160');

end