%noise sensitivity final
clear
rng(0)

eta_list = 0:0.05:1;
Eb_No_db_list = [11.5 12.5 13];


BER_eta = zeros(length(Eb_No_db_list), length(eta_list));
FER_eta = zeros(length(Eb_No_db_list), length(eta_list));
IT_eta  = zeros(length(Eb_No_db_list), length(eta_list));  % optional


T=60; %max decoder iteration

%eta=1; %.14
w=25; %.4; %25

flip_num=5;

p=4; % number of bits per symbol
q = 2^p; 

load('arith_16.mat');
load('204.102.3.6.16.mat'); 
load('G_204_102.mat');

h = full(h); % H matrix
N = size(h,2); %length of codeword
M = size(h,1); %number of parity checks
K=N-M; % msg length
R=K/N; % code rate

CN_lst=cell(M,1);

for i=1:M
    CN_lst{i,1}=find(h(i,:));
end

info_seq=randi([0 q-1], 1, K);

code_seq = gf_mat_mul(info_seq,G, add_mat, mul_mat); %given codeword
%Syndromes = decod_prod(code_seq,h,CN_lst, mul_mat, add_mat);

%initializing vector size
y=zeros(1,N);
hard_d_cmplx=zeros(1,N);
hard_d_gf16=zeros(1,N);


% corresponding QAM values as a complex array
qam16 = [
    -3 + 3i;  % 0011
    -3 + 1i;  % 0010
    -3 - 1i;  % 0001
    -3 - 3i;  % 0000
    -1 + 3i;  % 0111
    -1 + 1i;  % 0110
    -1 - 1i;  % 0101
    -1 - 3i;  % 0100
    1 + 3i;  % 1011
    1 + 1i;  % 1010
    1 - 1i;  % 1001
    1 - 3i;  % 1000
    3 + 3i;  % 1111
    3 + 1i;  % 1110
    3 - 1i;  % 1101
    3 - 3i   % 1100
    ];

qam_binary_map = [ 0  0  0  0;   % -3-3j
                   0  0  0  1;   % -3-1j
                   0  0  1  1;   % -3+1j
                   0  0  1  0;   % -3+3j
                   0  1  1  0;   % -1+3j
                   0  1  1  1;   % -1+1j
                   0  1  0  1;   % -1-1j
                   0  1  0  0;   % -1-3j
                   1  1  0  0;   % +1-3j
                   1  1  0  1;   % +1-1j
                   1  1  1  1;   % +1+1j
                   1  1  1  0;   % +1+3j
                   1  0  1  0;   % +3+3j
                   1  0  1  1;   % +3+1j
                   1  0  0  1;   % +3-1j
                   1  0  0  0];  % +3-3j


avg_pow = qam16'*qam16/q; %sum of squared magnitudes of all symbols(the total power)/q.
nrm_fct=sqrt(avg_pow); %normalization factor
gf16 = (0:q-1); %GF field symbols
alph_bin =  logical(fliplr(dec2bin(gf16, p) - 48)); % symbols in binary


targetFE = 1500;
max_gen  = 2e5;
                

for s = 1:length(Eb_No_db_list)

    EbN0 = Eb_No_db_list(s);
    Eb_No_linear = 10^(EbN0/10);

    

    for k = 1:length(eta_list)

        eta = eta_list(k);

        FE = 0;
        genFrame = 0;
        BE_Coded_total = 0;
        iterSum = 0;


         while(FE < targetFE && genFrame < max_gen)
               genFrame = genFrame + 1;
            
   
             c(1,1:N) = qam16(code_seq'+1,1);
    
                Eb_No_linear = 10^(EbN0/10);
                No = 1/(p * Eb_No_linear * R);
                sigma0 = sqrt(No/2)*nrm_fct;
                nse_std = eta * sigma0;
    
                n = sigma0*randn(1,N) + sigma0*randn(1,N)*1i;
                y = c + n;
    
            for j=1:N
                distance=abs(qam16-y(j));
                [~,min_idx]= min(distance);
                hard_d_cmplx(j) = qam16(min_idx); % hard decision in complex
                hard_d_gf16(j) = gf16(min_idx);
            end
           
            
            %calling the decoing function    
            [seqgf,failed,l]= decodeQamMinDis(code_seq,hard_d_cmplx, hard_d_gf16,...
            qam16, gf16,y, h, N, M, T, w,add_mat,mul_mat,div_mat,...
            CN_lst, nse_std,qam_binary_map,flip_num,No);  
                   
            iterSum = iterSum + l;
    
            %bit error calculation for decoded sequence 
        
            errors_coded_bit = zeros(1, K);
                for g = 1 : K            
                    if seqgf(g)~=code_seq(g)
                        s1 = dec2bin(code_seq(g),p);
                        s2 = dec2bin(seqgf(g),p);
                        code_seq_ = double(s1);
                        dec_seq_ = double(s2);
                        num_diff_bit = sum(code_seq_ ~= dec_seq_);
                        errors_coded_bit(g) = errors_coded_bit(g) + num_diff_bit;                
                    end     
                end
    
            bit_error = sum(errors_coded_bit);
            BE_Coded_total = BE_Coded_total + bit_error; % total no. of bit error in total frame generated in the current Eb/No
    
            %iters_cnr(i) = iters_cnr(i)+l;
    
            if (failed>0)
                FE=FE+1; 
            end       
         end


        BER_eta(s,k) = BE_Coded_total / (genFrame * N * p);
        FER_eta(s,k) = FE / genFrame;
        IT_eta(s,k)  = iterSum / genFrame;

        fprintf('Eb/No=%.1f dB | eta=%.2f | BER=%.3e | FER=%.3e | AvgIt=%.2f | Frames=%d\n', ...
            EbN0, eta, BER_eta(s,k), FER_eta(s,k), IT_eta(s,k), genFrame);

    end
end

figure; hold on; grid on;
set(gca,'YScale','log');
set(gca,'FontSize',12);      % tick numbers size
set(gcf,'Color','w');

markers = {'o','s','*'};
cols = lines(length(Eb_No_db_list));

for s = 1:length(Eb_No_db_list)
    plot(eta_list, BER_eta(s,:), ...
        markers{s}, ...
        'LineStyle','none', ...
        'MarkerSize',9, ...
        'LineWidth',2, ...
        'Color', cols(s,:), ...
        'MarkerFaceColor','none');
end

xlabel('Noise scaling factor \eta', 'FontSize', 18);
ylabel('BER', 'FontSize', 18);

legend(arrayfun(@(x) sprintf('E_b/N_0 = %.1f dB', x), Eb_No_db_list, ...
    'UniformOutput', false), ...
    'Location','best', 'FontSize',8);
