clear
rng(0)

Eb_No_db =10:1:15;

T=960; %max decoder iteration


eta=1; %.14 
w=25; %.4; %15

flip_num=5;

p=4; % number of bits per symbol
q = 2^p; 

load('arith_16.mat');
load('nbldpc_45_2_3_gf16_seed7.mat'); 


%h = full(h); % H matrix
N = size(h,2); %length of codeword
M = size(h,1); %number of parity checks
%K=N-M; % msg length
R=K/N; % code rate

CN_lst = str_cn_vn;


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

%initializing vector size

y = zeros(1, N);
hard_d_cmplx = zeros(1, N);
hard_d_gf16  = zeros(1, N);


FE=zeros(length(Eb_No_db),1);
genFrame=zeros(length(Eb_No_db),1);
iters_cnr=zeros(length(Eb_No_db),1);
BE_Coded=zeros(length(Eb_No_db),1);
BE_unCoded=zeros(length(Eb_No_db),1);


targetFE=500; %maximum FER to be observed
max_gen=1e6; % maximum number of frame to be generated

info_seq=randi([0 q-1], 1, K);

%code_seq = gf_mat_mul(info_seq,G, add_mat, mul_mat); %given codeword

code_seq = gf_vec_mat_mul(info_seq, G, add_mat, mul_mat);


Syndromes = decod_prod(code_seq, h, CN_lst, mul_mat, add_mat);

for i = 1:length(Eb_No_db)


    % Set the parameters based on the current Eb_N0 value
    %dynamic adjusting of number of flips to be done
    if Eb_No_db(i) >= 9 && Eb_No_db(i)<= 10
        flip_numd=5;      
    elseif Eb_No_db(i) > 10 && Eb_No_db(i)<= 11
        flip_numd=5;
    elseif Eb_No_db(i) > 11 && Eb_No_db(i)<= 14
        flip_num=5;
    end


     while(FE(i) < targetFE && genFrame(i)<max_gen)

        genFrame(i)=genFrame(i)+1;

        c(1,1:N) = qam16(code_seq'+1,1); % codeword in complex
        avg_symbol_energy = 1;

        Eb_No_linear(i)= 10.^(Eb_No_db(i)/10);
        No = avg_symbol_energy ./ (p * Eb_No_linear(i)*R); %noise spectral density
        sigma0 = sqrt(No/2)*nrm_fct ; %noise standard deviation
        nse_std=eta*sigma0; %noise perterbation used in flipping function E

        n = sigma0*randn(1,N)+sigma0*randn(1,N)*1i; %complex noise
        y = c + n; % codeword+noise, channel information

        for j=1:N
            distance=abs(qam16-y(j));
            [~,min_idx]= min(distance);
            hard_d_cmplx(j) = qam16(min_idx); % hard decision in complex
            hard_d_gf16(j) = gf16(min_idx);
        end
        
        errors = hard_d_cmplx ~= c; 
        n_errors_hard =sum(errors); %total number of symbol errors after hard decision made
              
        errors_uncoded_bit = zeros(1, K);

        %bit error calculation in hard decision
        for e = 1 : K          
                if hard_d_gf16(e)~=code_seq(e) 
                    s1 = dec2bin(code_seq(e),p);
                    s2 = dec2bin(hard_d_gf16(e),p);
                    code_seq_ = double(s1);
                    dec_seq_ = double(s2);
                    num_diff_bit = sum(code_seq_ ~= dec_seq_);
                    errors_uncoded_bit(e) = errors_uncoded_bit(e) + num_diff_bit;                
                end
        end
        
        un_bit_error = sum(errors_uncoded_bit); % no of bit error in each frame
        BE_unCoded(i)=BE_unCoded(i)+un_bit_error; % total no. of bit error in total frame generated in current Eb/No
        
        %calling the decoing function    
        [seqgf,failed,l]= decodeQamMinDis(code_seq,hard_d_cmplx, hard_d_gf16,...
        qam16, gf16,y, h, N, M, T, w,add_mat,mul_mat,div_mat,...
        CN_lst, nse_std,qam_binary_map,flip_num,No);  
               
        iters_cnr(i)=iters_cnr(i)+l;

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

        bit_error = sum(errors_coded_bit); % no. of bit error in each frame
        BE_Coded(i)=BE_Coded(i)+bit_error; % total no. of bit error in total frame generated in the current Eb/No

        %iters_cnr(i) = iters_cnr(i)+l;

        if (failed>0)
            FE(i)=FE(i)+1; 
        end       
     end
     fprintf('Eb/No = %.1f dB: BER (coded) = %.6e, BER (uncoded) = %.6e\n', ...
        Eb_No_db(i), ...
        BE_Coded(i) / (genFrame(i) * N * p), ...
        BE_unCoded(i) / (genFrame(i) * N * p));
end

BERunCoded= BE_unCoded ./(genFrame * N *p); %bit error rate for uncoded

BERCoded = BE_Coded ./ (genFrame * N *p); %bit error rate for coded


 figure;
 
 semilogy(Eb_No_db, BERCoded, 'gx-', 'LineWidth', 1.2);
 grid on;
 hold on;

 semilogy(Eb_No_db, BERunCoded, 'rx-','LineWidth', 1.2);
 hold off;

 ylim([10e-8 10e-1]);
 xlim([5 25]);
 xlabel('E_b/N_0 (dB)', 'FontSize', 14);
 ylabel('BER', 'FontSize', 14);
 %title('BER curve LDPC code (45,2,3) over GF(16)');
 hold off;
legend( 'NGDSF', 'Uncoded','Location', 'northeast');











