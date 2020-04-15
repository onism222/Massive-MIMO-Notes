%% Define parameters
%Empty workspace and close figures
close all;
clear;

%Define the range of BS antennas
Mvalues = [1, 10, 100];

%Angle of the desired UE
varphiDesired = pi/6;

%Range of angles of the interfering UE
varphiInterfererDegrees = [-180 : 1 : 180];
varphiInterfererRadians = varphiInterfererDegrees * (pi/180);

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Preallocate matrix for storing the simulation results
gfunction = zeros(length(varphiInterfererDegrees), length(Mvalues));


%% Go through all number of antennas
for m = 1 : length(Mvalues)
    
    %Generate channel response for the desired UE using (1.23)
    %assume \beta_0^0 = \beta_1^0 = 1
    hdesired = exp(1i * 2 * pi * antennaSpacing * sin(varphiDesired) * (0 : Mvalues(m) - 1)');
    
    %Go through all angles of interfering UE
    for n = 1 : length(varphiInterfererRadians)
        
        %Generate channel response for the interfering UE using (1.23)
        hinterfering = exp(1i * 2 * pi * antennaSpacing * sin(varphiInterfererRadians(n)) * (0 : Mvalues(m) - 1)');
        
        %Compute the g-function in (1.28), using its definition
        gfunction(n, m) = abs(hdesired' * hinterfering) ^ 2 / Mvalues(m);
        
    end
    
end

%% Plot the simulation results
figure;
hold on; box on;

plot(varphiInterfererDegrees, gfunction(:, 1), 'k-', 'LineWidth', 1);
plot(varphiInterfererDegrees, gfunction(:, 2), 'r--', 'LineWidth', 1);
plot(varphiInterfererDegrees, gfunction(:, 3), 'b-.', 'LineWidth', 1);

xlabel('Angle of interfering UE [degree]');
ylabel('$$g(\varphi_0^0, \varphi_1^0)$$', 'Interpreter', 'latex');

set(gca, 'Yscale', 'log');
xlim([-180, 180]);
ylim([1e-5, 1e2]);

legend('M=1', 'M=10', 'M=100', 'Location', 'NorthWest');
grid on
set(gca, 'color',  [1, 0.9, 0.8]);