% Eb/N0 points (dB)
Eb_No_db = [10 11 12 13 14 15];

% Coded BER values
BER_coded = [ ...
    7.525899e-03, ...
    3.085520e-03, ...
    6.687786e-04, ...
    3.151170e-04, ...
    8.095527e-05, ...
    7.478632e-06 ];

figure;
semilogy(Eb_No_db, BER_coded, 'o-', 'LineWidth', 1.8, 'MarkerSize', 7);
grid on;

xlabel('E_b/N_0 (dB)', 'FontSize', 14);
ylabel('Coded BER', 'FontSize', 14);
title('NB-LDPC (45,2,3) over 16-QAM — Coded BER', 'FontSize', 14);

xlim([9.5 15.5]);
ylim([1e-6 1e-2]);
