clc
close all
clear

%% Importing data

data = readtable("Airfoil_Analysis_Data - Spr26.xlsx");
atmosdata = [data(:,2) data(:,3)];

col1 = atmosdata{:,1};   % returns numeric array from first column
col2 = atmosdata{:,2};   % returns numeric array from second column
density = col1 ./ (col2 .* 287);

dynamicPressure = data(:,4);
freeV = sqrt(2 .* density .* data(:,4));

pressurecoef = ()