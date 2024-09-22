function Data = readLifeTables(filename, flag)
% Function to extract data from the life tables saved in a .xls file
%
% INPUTS: 
% filename: Name of the file
% flag:     Flag variable:
%           -if flag==1, consider males only
%           -if flag==2, consider both males and females
%
% OUTPUT:
% Data:     Struct with extracted data

% Extract data on the age:
Data.x = xlsread(filename,flag,'A3:A119');
% Extract data on survivors:
Data.lx = xlsread(filename,flag,'B3:B119');
% Extract data on deaths:
Data.dx = xlsread(filename,flag,'C3:C119');
% Extract data on death probability (per one thousand):
Data.qx = xlsread(filename,flag,'D3:D119');
% Extract data on years lived:
Data.Lx = xlsread(filename,flag,'E3:E119');
% Extract data on probabilities on prospective of survival:
Data.Px = xlsread(filename,flag,'F3:F119');
% Extract data on life expectancy:
Data.ex = xlsread(filename,flag,'G3:G119');

end