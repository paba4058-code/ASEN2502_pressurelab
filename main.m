clc
close all
clear

%% Importing data

data = readtable("Airfoil_Analysis_Data - Spr26.xlsx");
data2 = readtable("PortLocations.xlsx");
alpha = data{:,1};
portloc = data2{:,1};
deltax = data2{:,2};
deltay = data2{:,3};
atmosdata = [data(:,2) data(:,3)];
portdata = [data{:,6} data{:,7} data{:,8} data{:,9} data{:,10} data{:,11} data{:,12} data{:,13} data{:,14} data{:,5} data{:,15} data{:,16} data{:,17} data{:,18} data{:,19} data{:,20} data{:,21}];

%% Separate variables
P = atmosdata{:,1};   % Atmospheric Pressure
T = atmosdata{:,2};   % Atmospheric Temperature
testP = data{:,5};    % Test section pressure
dynamicPressure = data{:,4};
pressurecoeftemp = zeros(size(portdata,1), 17);
totalLift = zeros(24,1);
visc = 1.74e-5;
chordL = 0.0889;
Nupper = zeros(24,1);
Nlower = zeros(24,1);
Aupper = zeros(24,1);
Alower = zeros(24,1);

%% Calculations
T_K = T + 273.15;            % Convert to Kelvin
density = P ./ (287 .* T_K);

densityavg = mean(density);
freeV = sqrt(2 .* dynamicPressure ./ density);
freeVavg = mean(freeV);
ReNumber = densityavg * freeVavg * chordL / visc;

% Normal and axial trapezoidal approximations
numTests = size(portdata,1);

for j = 1:numTests

    Nupper(j) = 0;
    Aupper(j) = 0;
    Nlower(j) = 0;
    Alower(j) = 0;

    % Upper surface
    for i = 1:10
        p1 = portdata(j,i);
        p2 = portdata(j,i+1);

        Nupper(j) = Nupper(j) - 0.5*(p1 + p2)*deltax(i);
        Aupper(j) = Aupper(j) - 0.5*(p1 + p2)*deltay(i);
    end

    %  Lower surface
    for i = 10:16
        p1 = portdata(j,i);
        p2 = portdata(j,i+1);

        Nlower(j) = Nlower(j) - 0.5*(p1 + p2)*deltax(i);
        Alower(j) = Alower(j) - 0.5*(p1 + p2)*deltay(i);
    end
end


Ntotal = Nupper + Nlower;
Atotal = Aupper + Alower;

totalLift = Ntotal .* cosd(alpha) - Atotal .* sind(alpha);

% coeff of pressure calculation

for i=1:17
pressurecoeftemp(:,i) = (portdata(:,i) - testP) ./ dynamicPressure;
end

pressurecoeftemp(:,10) = 0;
pressurecoef = pressurecoeftemp';

% coeff of lift calc

liftcoeff1 =(dynamicPressure * chordL);
liftcoeff = totalLift ./ liftcoeff1;

%% Post Processing
numTests = size(pressurecoef, 2);

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

    upperCp = pressurecoef(upperIdx, k);
    lowerCp = pressurecoef(lowerIdx, k);

    upperX = portloc(upperIdx);
    lowerX = portloc(lowerIdx);

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

figure
plot(alpha, liftcoeff, 'o-','LineWidth',1.5)
xlabel('\alpha (deg)')
ylabel('Lift (pressure integrated)')
grid on


