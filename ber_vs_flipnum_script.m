clear; rng(0);

flip_list = [3 4 5 6 7];
Eb_No_db_list = [10.5 11.5 12.5];

w = 25;
eta = 1;

BER_f = zeros(length(Eb_No_db_list), length(flip_list));
FER_f = zeros(length(Eb_No_db_list), length(flip_list));
IT_f  = zeros(length(Eb_No_db_list), length(flip_list));

T = 500;

p=4; q=2^p;

load('arith_16.mat');
load('204.102.3.6.16.mat'); 
load('G_204_102.mat');

h = full(h);
N = size(h,2);
M = size(h,1);
K=N-M;
R=K/N;

y = zeros(1,N);
hard_d_cmplx = zeros(1,N);
hard_d_gf16  = zeros(1,N);

CN_lst=cell(M,1);
for i=1:M
    CN_lst{i,1}=find(h(i,:));
end

qam16 = [
    -3+3i;-3+1i;-3-1i;-3-3i;
    -1+3i;-1+1i;-1-1i;-1-3i;
     1+3i; 1+1i; 1-1i; 1-3i;
     3+3i; 3+1i; 3-1i; 3-3i];

qam_binary_map = [
0 0 0 0; 0 0 0 1; 0 0 1 1; 0 0 1 0;
0 1 1 0; 0 1 1 1; 0 1 0 1; 0 1 0 0;
1 1 0 0; 1 1 0 1; 1 1 1 1; 1 1 1 0;
1 0 1 0; 1 0 1 1; 1 0 0 1; 1 0 0 0];

avg_pow = qam16'*qam16/q;
nrm_fct=sqrt(avg_pow);
gf16 = (0:q-1);

targetFE = 300;
max_gen  = 2e5;

info_seq = randi([0 q-1], 1, K);
code_seq = gf_mat_mul(info_seq, G, add_mat, mul_mat);

for s = 1:length(Eb_No_db_list)

    EbN0 = Eb_No_db_list(s);
    Eb_No_linear = 10^(EbN0/10);

    for k = 1:length(flip_list)

        flip_num = flip_list(k);

        FE = 0; genFrame = 0;
        BE_Coded_total = 0;
        iterSum = 0;

        

        while(FE < targetFE && genFrame < max_gen)

            genFrame = genFrame + 1;

            c(1,1:N) = qam16(code_seq'+1,1);

            No = 1/(p * Eb_No_linear * R);
            sigma0 = sqrt(No/2)*nrm_fct;
            nse_std = eta * sigma0;

            n = sigma0*randn(1,N) + sigma0*randn(1,N)*1i;
            y = c + n;

            for j = 1:N
                [~,min_idx] = min(abs(qam16 - y(j)));
                hard_d_cmplx(j) = qam16(min_idx);
                hard_d_gf16(j)  = gf16(min_idx);
            end

            [seqgf, failed, l] = decodeQamMinDis(code_seq, hard_d_cmplx, hard_d_gf16, ...
                qam16, gf16, y, h, N, M, T, w, add_mat, mul_mat, div_mat, ...
                CN_lst, nse_std, qam_binary_map, flip_num, No);

            iterSum = iterSum + l;

            % bit error
            for g = 1:K
                if seqgf(g) ~= code_seq(g)
                    s1 = dec2bin(code_seq(g),p);
                    s2 = dec2bin(seqgf(g),p);
                    BE_Coded_total = BE_Coded_total + sum(double(s1) ~= double(s2));
                end
            end

            if failed>0
                FE = FE + 1;
            end
        end

        BER_f(s,k) = BE_Coded_total / (genFrame * N * p);
        FER_f(s,k) = FE / genFrame;
        IT_f(s,k)  = iterSum / genFrame;

        fprintf('Eb/No=%.1f | flip=%d | BER=%.3e\n', EbN0, flip_num, BER_f(s,k));
    end
end

figure; hold on; grid on;
set(gca,'YScale','log');
set(gca,'FontSize',16);
set(gcf,'Color','w');

plot(flip_list, BER_f(1,:), 'o','MarkerSize',8,'LineWidth',1.5);
plot(flip_list, BER_f(2,:), 's','MarkerSize',8,'LineWidth',1.5);
plot(flip_list, BER_f(3,:), '*','MarkerSize',9,'LineWidth',1.5);

xlabel('flip\_num', 'FontSize', 22);
ylabel('BER', 'FontSize', 22);

legend('E_b/N_0 = 10.5 dB', ...
       'E_b/N_0 = 11.5 dB', ...
       'E_b/N_0 = 12.5 dB', ...
       'Location','best','FontSize',16);

print(gcf,'BER_vs_flipnum','-dpdf','-painters');

