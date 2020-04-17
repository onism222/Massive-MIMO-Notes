%% Define parameters
%Empty workspace and close figures
close all;
clear;

%Define the range of number of BS antennas
M = [10, 100];

%Select the number of Monte Carlo realizations of the Rayleigh fading
numberOfRealizations = 1000000;

%Generate random UE angles from 0 to 2*pi
varphiDesired = 2 * pi * rand(1, numberOfRealizations);
varphiInterfering = 2 * pi * rand(1, numberOfRealizations);

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Preallocate matrix for storing the simulation results
interferenceGainLoS = zeros(numberOfRealizations, length(M));

%Compute the argument, that appears in (1.28), for different UE angles
argument = (pi * antennaSpacing * (sin(varphiDesired(:)) - sin(varphiInterfering(:))));

%% Go through different number of antennas
for mindex = 1 : length(M)

    %Compute the g-function in (1.28)
    interferenceGainLoS(:, mindex) =  (sin(argument * M(mindex))).^2 ./ (sin(argument)).^2 / M(mindex);
    
end

CDFvalues_LoS = linspace(0, 1, numberOfRealizations); 

%Compute the CDF of the relative interference gain for NLoS using the
%Exp(1)-distribution
interferenceGainNLoS = logspace(-6, 2, 10000);
CDFvalues_NLoS = 1 - exp(-interferenceGainNLoS);

%% Plot the simulation results
figure;
hold on; box on;

plot(interferenceGainNLoS, CDFvalues_NLoS, 'r--', 'LineWidth', 1.5);
plot(sort(interferenceGainLoS(:, 1)), CDFvalues_LoS, 'b-.', 'LineWidth', 1.5);
plot(sort(interferenceGainLoS(:, 2)), CDFvalues_LoS, 'k-', 'LineWidth', 1.5);

set(gca, 'Xscale', 'log');
xlim([1e-6, 1e2]),
xlabel('Relative interference gain');
ylabel('CDF');

legend('NLoS, any M', 'LoS, M=10', 'LoS, M=100', 'Location', 'NorthWest');
grid on
set(gca, 'color',  [1, 0.9, 0.8]);
