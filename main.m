clc
close all
clear

%% Importing data

data = readtable("Airfoil_Analysis_Data - Spr26.xlsx");
data2 = readtable("PortLocations.xlsx");
portloc = data2{:,1};
atmosdata = [data(:,2) data(:,3)];
portdata = [data{:,6} data{:,7} data{:,8} data{:,9} data{:,10} data{:,11} data{:,12} data{:,13} data{:,14} data{:,5} data{:,15} data{:,16} data{:,17} data{:,18} data{:,19} data{:,20} data{:,21}];

%% Separate variables
P = atmosdata{:,1};   % Atmospheric Pressure
T = atmosdata{:,2};   % Atmospheric Temperature
testP = data{:,5};    % Test section pressure
dynamicPressure = data{:,4};
pressurecoeftemp = zeros(size(portdata,1), 17);

%% Calculations
T_K = T + 273.15;            % Convert to Kelvin
density = P ./ (287 .* T_K);
freeV = sqrt(2 .* dynamicPressure ./ density);

for i=1:17
pressurecoeftemp(:,i) = (portdata(:,i) - testP) ./ dynamicPressure;
end

pressurecoeftemp(:,10) = 0;
pressurecoef = pressurecoeftemp';

%% Post Processing
% Number of test cases (angles of attack)
numTests = size(pressurecoef, 2);

% Choose layout (roughly square)
nCols = ceil(sqrt(numTests));
nRows = ceil(numTests / nCols);

figure
tiledlayout(nRows, nCols, ...
    'TileSpacing','compact', ...
    'Padding','compact')
cpMin = min(pressurecoef,[],'all');
cpMax = max(pressurecoef,[],'all');
upperIdx = 1:10;      % upper surface ports
lowerIdx = 10:17;   % lower surface ports


for k = 1:numTests
    nexttile
    hold on

    h1 = plot(portloc(upperIdx), pressurecoef(upperIdx,k), 'r-o','LineWidth',1.2);
    h2 = plot(portloc(lowerIdx), pressurecoef(lowerIdx,k), 'b-s','LineWidth',1.2);
    set(gca,'YDir','reverse')
    ylim([cpMin cpMax])
    grid on

    title(sprintf('\\alpha = %g^\\circ', data{k,1}))  
    xlabel('x/c')
    ylabel('C_p')
    
    hold off
end

legend([h1 h2], 'Upper surface', 'Lower surface', 'Location','bestoutside')
