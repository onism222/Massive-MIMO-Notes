%Empty workspace and close figures
close all;
clear;


%Define the SNR range for analytical curves
SNRdB = -10 : 0.1 : 30;
SNR = 10 .^ (SNRdB/10);

%Define the SNR range for Monte Carlo simulations
SNRdB_montecarlo = -10 : 5 : 30;
SNR_montecarlo = 10 .^ (SNRdB_montecarlo/10); % SNR_montecarlo = \frac{p}{\sigma^2}

%Define the different beta_bar values (strength of inter-cell interference)
betabar = [1e-1, 1e-3]';

%Preallocate matrices for storing the simulation results
SE_LoS = zeros(length(betabar), length(SNR));
SE_NLoS = zeros(length(betabar), length(SNR));
SE_NLoS_montecarlo = zeros(length(betabar), length(SNR_montecarlo));

%Select number of Monte Carlo realizations of the Rayleigh fading
numberOfFadingRealizations = 100000;


%% Go through different strengths of the interference
for b = 1 : length(betabar)
    
    %Compute SE under line-of-sight (LoS) propagation as in (1.17)
    SE_LoS(b, :) = log2(1 + 1 ./ (betabar(b) + 1 ./ SNR)); 
    
    
    %Generate uncorrelated Rayleigh fading channel realizations
    fadingRealizationsDesired = (randn(numberOfFadingRealizations,1) + 1i * randn(numberOfFadingRealizations, 1)) / sqrt(2);%h_0^0
    fadingRealizationsInterference = (randn(numberOfFadingRealizations, 1) +1i * randn(numberOfFadingRealizations, 1)) / sqrt(2);%h_1^0
    
    %Compute SE under non-line-of-sight (NLoS) propagation from the first
    %line in (1.18), using Monte Carlo simulations for the channel realizations
    SE_NLoS_montecarlo(b, :) = mean(log2(1 + abs(fadingRealizationsDesired) .^ 2 * ...
        SNR_montecarlo ./ (abs(fadingRealizationsInterference) .^ 2 * SNR_montecarlo*betabar(b) +1)),1);
    %note that SNR_montecarlo = \frac{p}{\sigma^2}
    %M = mean(A,dim) 返回维度 dim 上的均值。例如，如果 A 为矩阵，则 mean(A,2) 是包含每一行均值的列向量。
    
    
    %Compute SE under non-line-of-sight (NLoS) propagation as in (1.18)
    SE_NLoS(b,:) = (exp(1 ./ SNR) .* expint(1 ./ SNR) - exp(1 ./ (betabar(b) * SNR)) .* expint(1 ./ (betabar(b) * SNR)))/((1 - betabar(b)) * log(2));
    
end


%% Plot the simulation results
figure;
hold on; box on;

for b = 1 : length(betabar)
    
    plot(SNRdB, SE_LoS(b, :), 'k-', 'LineWidth', 1);
    plot(SNRdB_montecarlo, SE_NLoS_montecarlo(b, :), 'rd', 'LineWidth', 1);
    plot(SNRdB, SE_NLoS(b, :), 'b-.', 'LineWidth', 1);
    
end

xlabel('SNR [dB]');
ylabel('Average SE [bit/s/Hz]');

legend('beta = 0.1 and 0.001, Analytical LoS', 'beta = 0.1 and 0.001, Montecarlo NLoS',...
    'beta = 0.1 and 0.001, Analytical NLoS', 'Location','NorthWest');
ylim([0 10]);
grid on
set(gca, 'color',  [1, 0.9, 0.8]);
