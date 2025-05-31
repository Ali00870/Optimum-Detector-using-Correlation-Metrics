
 clc; clear; close all;
function s = pam_gray_map(n, M, Delta)
%PAM_GRAY_MAP  Gray-coded M-PAM mapper
%   n     : vector of integers in {0,1,…,M-1}
%   M     : number of levels (must be a power Esof 2)
%   Delta : spacing between adjacent constellation points
%
%   s     : vector of mapped PAM symbols

    %--- 1. Compute Gray code of each integer: g = n XOR floor(n/2)
    g = bitxor(uint32(n), bitshift(uint32(n), -1));

    %--- 2. Map Gray index to amplitude (centered at zero)
    %    levels: {-(M-1)*Delta/2, …, +(M-1)*Delta/2}
    s = Delta * (double(g) - (M-1)/2);

end
% Parameters
M = 4;
k = log2(M);
Delta = 2;
num_symbols = 1e5;
T = 1;
dt = 1 / (1);
samples_per_symbol = round(T / dt);
EbN0_dB = -10:1:10;
BER_simulated = zeros(size(EbN0_dB));

% PAM symbol mapping
waveforms = 0:M-1;
gray_waveforms = pam_gray_map(waveforms, M, Delta);
[gray_waveforms, gray_idx] = sort(gray_waveforms);
bitmapping = waveforms(gray_idx);

% Average symbol energy
Es_per_level = gray_waveforms.^2 * T;
Es_avg = mean(Es_per_level);
Eb_avg = Es_avg / k;

% Generate random symbols
n = randi([0 M-1], 1, num_symbols);
transmitted_amps = pam_gray_map(n, M, Delta);
transmitted_symbols = repmat(transmitted_amps', 1, samples_per_symbol);
EbN0_linear = 10.^(EbN0_dB/10);
% Theoretical BER (for M-PAM with Gray coding in AWGN)
theoretical_BER = berawgn(EbN0_dB, 'pam', M, 'binary');
for idx = 1:length(EbN0_dB)
    % Noise variance for this Eb/N0
    N0 = Eb_avg / (10^(EbN0_dB(idx)/10))/dt;

    % Generate AWGN
    noise = sqrt(N0/2) * randn(num_symbols, samples_per_symbol);

    % Received signals
    received_symbols = transmitted_symbols + noise;

    num_bit_errors = 0;

    for i = 1:num_symbols
        received = received_symbols(i, :);
        metrics = zeros(1, M);

        % Correlation receiver decision
        for j = 1:M
            template = gray_waveforms(j) * ones(1, samples_per_symbol);
            correlation = 2 * sum(received .* template) * dt;
            energy = gray_waveforms(j)^2 * T;
            metrics(j) = correlation - energy;
        end

        [~, detected_idx] = max(metrics);
        detected_symbol_index=detected_idx;
        true_symbol_index=find(gray_waveforms== transmitted_symbols(i,1));

       detected_gray=bitmapping(detected_symbol_index);
       true_gray=bitmapping(true_symbol_index);

        % Bits
        true_bits = de2bi(true_gray, k, 'left-msb');
        detected_bits = de2bi(detected_gray, k, 'left-msb');

        % Count bit errors
        num_bit_errors = num_bit_errors + sum(true_bits ~= detected_bits);
    end

    % BER for this Eb/N0
    BER_simulated(idx) =num_bit_errors / (num_symbols * k);
    disp(['Eb/N0 = ', num2str(EbN0_dB(idx)), ' dB, Simulated BER = ', num2str(BER_simulated(idx)),' Theoretical BER= ',num2str(theoretical_BER(idx))])
end

% Plot results
figure;
semilogy(EbN0_dB, BER_simulated, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor','b');
hold on;
semilogy(EbN0_dB, theoretical_BER, 'r--', 'LineWidth', 2);
grid on;
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate');
title('BER vs Eb/N0 for Gray-coded 4-PAM with Matched Filter Receiver');
legend('Simulated', 'Theoretical', 'Location', 'southwest');
axis([min(EbN0_dB) max(EbN0_dB) 1e-8 1]);
