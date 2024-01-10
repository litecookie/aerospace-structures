% Matlab code for Aerospace Structures Project - 1

% ------------------------------------------------------------------------
% Part A - Failure Mechanism Map

% Define tbyc values from 0.01 to 0.30 with a 0.01 increment
tbyc = 0.01:0.01:0.30;

% Initialize arrays for failure mode intersections
mccbyl = zeros(1, 30);
cicbyl = zeros(1, 30);
micbyl = zeros(1, 30);

% Define face sheet properties
rhof = 1600 + (30 * 1);
Ef = (40 + 88) * (10^9);
sigmaf = (200 + (100 * 9)) * (10^6);

% Define foam core properties
rhoc = 20 + (5 * 67);
sigmac = (0.5 + 6.5 * ((67 / 100)^(1.5))) * (10^6);
tauc = (0.5 + 4.5 * ((67 / 100)^(1.5))) * (10^6);

% Calculate constants to use in the equations
cnum = pi^2 * Ef * sigmac^2;
cden = 192 * sigmaf^3;
cden2 = 24 * tauc^3;

K = cnum / cden;
K2 = cnum / cden2;

for i = 1:30
    % Microbuckling and Core Shear
    mccbyl(i) = tauc / (2 * sigmaf * tbyc(i));
    
    % Microbuckling and Indentation
    tbycc = (tbyc(i) + 1)^2;
    micbyl(i) = sqrt(K / tbycc);
    
    % Core Shear and Indentation
    cicbyl(i) = tbycc / (K2 * tbyc(i)^3);
end

% Create a meshgrid
num_values = 1000;
newtbyc = linspace(0, 0.30, num_values);
newcbyl = linspace(0, 0.30, num_values);
[ny, nx] = meshgrid(newtbyc, newcbyl);

% Calculate Pmcap, Pccap, and Picap
Pmcap = 4 .* ny .* nx .* (ny + 1) .* nx;
Pccap = 2 .* (tauc / sigmaf) .* (ny + 1) .* nx;
Ki = ((pi^2 .* Ef .* sigmac^2 / 3).^(1/3)) ./ sigmaf;
Picap = Ki .* ny .* nx .* (ny + 1).^(1/3) .* nx.^(1/3);

% Find the minimum values and their corresponding indices
min_values = min(min(Pmcap, Pccap), Picap);
index_matrix = zeros(size(min_values));

for i = 1:num_values
    for j = 1:num_values
        if min_values(i, j) == Pmcap(i, j)
            index_matrix(i, j) = 1;
        elseif min_values(i, j) == Pccap(i, j)
            index_matrix(i, j) = 2;
        else
            index_matrix(i, j) = 3;
        end
    end
end

% ------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Part b - Optimal Design Trajectory

%Failure due to microbuckling
rhobar = rhoc/rhof;
ymb = rhobar / (2 * (1 - rhobar));

%Failure due to elastic indentation - the line plotted is out of bounds.
yei = ( 3 * rhobar ) / (2 - 4 * rhobar);

% Point of intersection between microbuckling optimal design curve and
% microbuckling - indentation failure plot pair.
tbycc = (ymb + 1)^2;
intersection_x1 = sqrt(K / tbycc);

% ------------------------------------------------------------------------

% -------------------------------------------------------------------------
% PART C - Minimal Mass vs Load 

Pcap = 0.00001:0.00001:0.001;

Ebar = Ef/sigmaf;
sigmabar = sigmac/sigmaf;

den = 9*(pi^2)*(sigmabar^2)*(Ebar);
rhoe = rhobar*((2 - rhobar)^3);

Kei = 4 * ((rhoe/den)^(1/4));

Mcapei = Kei.*(Pcap.^(3/4));
Mcapmb = sqrt(Pcap.*(rhobar*(2 - rhobar)));

% Find the point of intersection
intersection_x2 = Pcap(abs(Mcapei - Mcapmb) < 0.0001);  % Adjust the tolerance as needed

Pcap1 = Pcap(Pcap <= intersection_x2);  % x-values for the first segment
Pcap2 = Pcap(Pcap >= intersection_x2);  % x-values for the second segment

% corresponding y-values 
Mcapei1 = Kei.*(Pcap1.^(3/4));
Mcapmb1 = sqrt(Pcap2.*(rhobar*(2 - rhobar)));

% ------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Part D - Example Beam Optimal Design

Pgiven = 75 * 10^3;
Lgiven = 2.5;
bgiven = 200 * 10^-3;

Pf = Pgiven / (bgiven * Lgiven * sigmaf);

% The value shows that for the load condition, minimum mass condition falls
% in the microbuckling region. (~0.13 * 10^-3)

Mf = sqrt(Pf*(rhobar*(2 - rhobar)));

cf = 0.01:0.01:0.08;

tf = ( Lgiven * Mf - rhobar .* cf ) / 2;

% ------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Plot A 1 - Failure Region Plot

figure(1);
hold on; 

plot(mccbyl, tbyc, 'b');
plot(micbyl, tbyc, 'g');
plot(cicbyl(5:end), tbyc(5:end), 'r');

% Additional plot settings
xlabel('C/L');
ylabel('T/C');
title('Project 1 Plot A 1');
legend('Core Shear - Elastic Indentation','Microbuckling - Elastic Indentation', 'Core Shear - Elastic Indentation');

% ------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Plot A 2 - Contour Plot

figure(2);
hold on;

plot(mccbyl, tbyc, 'b');
plot(micbyl, tbyc, 'g');
plot(cicbyl(5:end), tbyc(5:end), 'r');

custom_colormap = [0 1 0; 0 0 1; 1 0 0]; 
colormap(custom_colormap);  % Apply the custom colormap

contourf(nx, ny, index_matrix, 'LineStyle', 'none');

% Additional plot settings
xlabel('C/L');
ylabel('T/C');
title('Project 1 Plot A 2');
legend('Core Shear','Microbuckling', 'Elastic Indentation');

% ------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Plot B - Optimal Design Plots

figure(3);
hold on;

plot(mccbyl, tbyc, 'b');
plot(micbyl, tbyc, 'g');
plot(cicbyl(5:end), tbyc(5:end), 'r');

custom_colormap = [0 1 0; 0 0 1; 1 0 0]; 
colormap(custom_colormap);  % Apply the custom colormap

contourf(nx, ny, index_matrix, 'LineStyle', 'none');

% Plot the first part of optimal design curve 
% plot([0, intersection_x1], [ymb, ymb], 'black', 'LineWidth', 2);  

% Plot for both the lagrange multiplier calculated design trajectories
plot([0, 0.3], [ymb, ymb], 'yellow', 'LineWidth', 2);  
plot([0, 0.3], [yei, yei], 'black', 'LineWidth', 2);  

% Additional plot settings
xlabel('C/L');
ylabel('T/C');
title('Project 1 Plot B');
legend('Core Shear','Microbuckling', 'Elastic Indentation');

% ------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Plot C - Mass Minima Plots

figure(4);
hold on;
plot(Pcap, Mcapmb, 'r', 'LineWidth',1);
plot(Pcap, Mcapei, 'b','LineWidth',1);

% Plot the curve segments
plot(Pcap1, Mcapei1, 'blue', 'LineWidth', 3);  % Plot the first segment of curve 1 in blue
hold on;  % Keep the current figure active
plot(Pcap2, Mcapmb1, 'red', 'LineWidth', 3);  % Plot the second segment of curve 2 in red


% Additional plot settings
xlabel('Pcap');
ylabel('Mcap');
title('Project 1 Plot C');
legend('Microbuckling', 'Elastic Indentation');

% ------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Part D - Example Beam Geometry Plot

figure(5);
hold on;
plot(cf, tf, 'r');

% Additional plot settings
xlabel('C(in m)');
ylabel('T(in m)');
title('Project 1 Plot D');

hold off;
