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


deltax = deltax / chordL;
deltay = deltay / chordL;


%% Calculations
T_K = T + 273.15;
density = P ./ (287 .* T_K);

densityavg = mean(density);
freeV = sqrt(2 * dynamicPressure ./ density);
freeVavg = mean(freeV);
ReNumber = densityavg * freeVavg * chordL / visc;

% coeff of pressure calculation

for i=1:17
pressurecoeftemp(:,i) = (portdata(:,i) - testP) ./ dynamicPressure;
end

pressurecoeftemp(:,10) = 0;
pressurecoef = pressurecoeftemp';


% Normal and axial trapezoidal approximations
numTests = size(portdata,1);

for j = 1:numTests

    Nupper(j) = 0;
    Aupper(j) = 0;
    Nlower(j) = 0;
    Alower(j) = 0;

    % Upper surface
    for i = 1:10
        p1 = pressurecoef(i,j);
        p2 = pressurecoef(i+1,j);

        Nupper(j) = Nupper(j) - 0.5*(p1 + p2)*deltax(i);
        Aupper(j) = Aupper(j) + 0.5*(p1 + p2)*deltay(i);
    end

    %  Lower surface
    for i = 10:17-1

        p1 = pressurecoef(i,j);
        p2 = pressurecoef(i+1,j);

        Nlower(j) = Nlower(j) + 0.5*(p1 + p2)*deltax(i);
        Alower(j) = Alower(j) - 0.5*(p1 + p2)*deltay(i);
    end

    Nupper(j) = Nupper(j) * dynamicPressure(j);
Aupper(j) = Aupper(j) * dynamicPressure(j);
Nlower(j) = Nlower(j) * dynamicPressure(j);
Alower(j) = Alower(j) * dynamicPressure(j);

end


Ntotal = Nupper + Nlower;
Atotal = Aupper + Alower;

totalLift = Ntotal .* cosd(alpha) - Atotal .* sind(alpha);


% coeff of lift calc

liftcoeff1 =(dynamicPressure * chordL);
liftcoeff = totalLift ./ liftcoeff1;

%% Post Processing
numTests = 3;

figure
tiledlayout(1, numTests, ...
    'TileSpacing','compact', ...
    'Padding','compact')

cpMin = min(pressurecoef,[],'all');
cpMax = max(pressurecoef,[],'all');

% Wrapped Cp and x/c (all tests at once)
cpWrap = [
    pressurecoef(1:9, :);
    pressurecoef(10, :);
    pressurecoef(11:17, :);
    pressurecoef(1, :)
];

xWrap = [
    portloc(1:9);
    portloc(10);
    portloc(11:17);
    portloc(1)
];

% alpha ≈ -3 deg (near zero lift)
k = 6;
nexttile
hold on
plot(xWrap(1:10,:), cpWrap(1:10,k), 'r-o','LineWidth',1.2)
plot(xWrap(10:18,:), cpWrap(10:18,k), 'b-s','LineWidth',1.2)
set(gca,'YDir','reverse')
ylim([cpMin cpMax])
grid on
title('Angle of Attack = -3 (Lift \approx 0)')
xlabel('x/c')
ylabel('C_p')

% alpha = 8 deg
k = 17;
nexttile
hold on
plot(xWrap(1:10,:), cpWrap(1:10,k), 'r-o','LineWidth',1.2)
plot(xWrap(10:18,:), cpWrap(10:18,k), 'b-s','LineWidth',1.2)
set(gca,'YDir','reverse')
ylim([cpMin cpMax])
grid on
title('Angle of Attack = 8')
xlabel('x/c')
ylabel('C_p')

% alpha = 10 deg (post-stall)
k = 19;
nexttile
hold on
plot(xWrap(1:10,:), cpWrap(1:10,k), 'r-o','LineWidth',1.2)
plot(xWrap(10:18,:), cpWrap(10:18,k), 'b-s','LineWidth',1.2)
set(gca,'YDir','reverse')
ylim([cpMin cpMax])
grid on
title('Angle of Attack = 10 (Post-stall)')
xlabel('x/c')
ylabel('C_p')
legend('Upper Surface','Lower Surface')

figure
plot(alpha, liftcoeff, 'o-','LineWidth',1.5)
title('Coefficient of Lift over 24 angle of attack tests')
xlabel('\alpha (deg)')
ylabel('Lift coefficient')
grid on

figure
plot(alpha, totalLift, 'o-','LineWidth',1.5)
xlabel('\alpha (deg)')
ylabel('Lift')
grid on
