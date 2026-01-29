clear; clc;

load('arith_16.mat');   % add_mat, mul_mat, div_mat

% Generate H
[h, str_cn_vn, Hb] = gen_nbldpc_45_2_3_for_decod_prod(7);
disp('H generated'); disp(size(h));

% Generate G
[G, K] = gf16_generator_from_H(h, add_mat, mul_mat, div_mat);
fprintf('Code dimension K = %d\n', K);

% Encode random message
info_seq = uint8(randi([0 15], 1, K));
code_seq = gf_vec_mat_mul(info_seq, G, add_mat, mul_mat);   % recommended

% Syndrome check
S = decod_prod(code_seq, h, str_cn_vn, mul_mat, add_mat);
disp('Syndrome (should be all zero):');
disp(S);
disp(['All zero = ', num2str(all(S==0))]);

% Save
save('nbldpc_45_2_3_gf16_seed7.mat','h','str_cn_vn','G','K','Hb');
