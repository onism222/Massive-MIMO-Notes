%Empty workspace and close figures
close all;
clear;

%Define the SNR range for analytical curves
SNRdB = -10 : 0.1 : 50;
SNR = 10 .^ (SNRdB/10);

%Define the different beta_bar values (strength of inter-cell interference)
betabar = [1e-1, 1e-3]';

%Preallocate matrices for storing the simulation results
SE_LoS = zeros(length(betabar), length(SNR));

%% Go through different strengths of the interference
for b = 1 : length(betabar)  
    
    %Compute SE under line-of-sight (LoS) propagation as in (1.17)
    SE_LoS(b, :) = log2(1 + 1 ./ (betabar(b) + 1 ./ SNR));  
    
end

%% Plot the simulation results
figure;
hold on; box on; 
   
plot(SNRdB, SE_LoS(1, :), 'r--', 'LineWidth', 1.5);
plot(SNRdB, SE_LoS(2, :), 'b-.', 'LineWidth', 1.5);

xlabel('SNR [dB]');
ylabel('Average SE [bit/s/Hz]');

legend('beta = 0.1, LoS', 'beta = 0.001, LoS', 'Location', 'SouthEast');
ylim([0, 10]);
grid on
set(gca, 'color',  [1, 0.9, 0.8]);
