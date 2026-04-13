% Eb/N0 values (dB)
EbN0 = [9 10 11 12 13];

% BER data from table (eta = 1, flip_num = 5)
BER_T60  = [1.476225e-02 7.522059e-03 2.305952e-03 5.055730e-04 1.449307e-04];
BER_T460 = [1.449265e-02 6.518825e-03 1.682763e-03 2.208417e-04 3.348334e-05];
BER_T600 = [1.471814e-02 6.688725e-03 1.790067e-03 1.797616e-04 7.382471e-06 ];
BER_T960 = [1.442157e-02 6.365196e-03 1.749088e-03 1.494500e-04 2.709264e-06 ];

% Create figure
figure;
semilogy(EbN0, BER_T60,  '-o', 'LineWidth', 1); hold on;
semilogy(EbN0, BER_T460, '-s', 'LineWidth', 1);
semilogy(EbN0, BER_T600, '-^', 'LineWidth', 1);
semilogy(EbN0, BER_T960, '-d', 'LineWidth', 1);

% Grid and labels
grid on;
xlabel('$E_b/N_0$ (dB)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);

% Title
title('BER Performance of NGDSF (GF16), w = 25, \eta = 1, flip\_num = 5', ...
      'FontSize', 12);

% Legend
legend('T = 60', 'T = 460', 'T = 600', 'T = 960', ...
       'Location', 'southwest');

% Axis limits for clean comparison
ylim([1e-7 1e-1]);
xlim([9 14]);

% Improve aesthetics
set(gca, 'FontSize', 11, 'LineWidth', 1);
