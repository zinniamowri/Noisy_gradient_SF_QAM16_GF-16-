function [seqgf,failed_init,l]= decodeQamMinDis(code_seq, hard_d_cmplx, hard_d_gf16,...
    qam16, gf16,y, h, N, M, T, w, add_mat,mul_mat,div_mat,...
    CN_lst, nse_std,qam_binary_map,flip_num, No)

	d_deco_0 = hard_d_cmplx;
    seqgf= hard_d_gf16;
    d_dec_f =hard_d_cmplx;
	l = 0;
    %nse_std2 = (nse_std^2)/2; %used as scaling factor
    
    d_gf16 = hard_d_gf16;

	S = decod_prod(hard_d_gf16, h, CN_lst, mul_mat, add_mat);
    Sb = double(S==0); % satisfied syndrome S=0 is now considered as 1
	num_satisfied = sum(Sb);
    failed_init= M-num_satisfied;
    % Initialize tracking variables for best codeword
    min_failed = failed_init;
 
    best_d_gf16 = d_gf16;
    
    Hb=double(h>0);

    while (l < T)
          l = l + 1;
          %disp(l);   
          Sb = double(S==0); % satisfied syndrome S=0 is now considered as 1
          WSH = w*Sb*Hb;

          %No=2*sigma^2, total noise variance in QAM16
          %(-abs(y - d_dec_f).^2/2*sigma^2), this part comes from the log
          %of PDF of QAM16, which relates to LLR equation

          E = (-abs(y - d_dec_f).^2/No) + WSH ; %+ nse_std*randn(1,N); 

          % Select most unreliable symbols, by given number flip_num
          [~, idx] = mink(E, flip_num);

          temp_gf16 = d_gf16;

          for j=1:length(idx)
            current_symbol = d_dec_f(idx(j));
            distance=abs(qam16-current_symbol);
            [~, sorted_idx]=mink(distance,6);
            closest_symbols=gf16(sorted_idx);
            rand_idx=randi([1,4]);
            new_symbol=closest_symbols(rand_idx);
            temp_gf16(idx(j))= new_symbol; 
          end

            S_temp=decod_prod(temp_gf16, h, CN_lst, mul_mat, add_mat);
            Sb_temp = double(S_temp==0); % satisfied syndrome S=0 is 1
	        num_satisfied_temp = sum(Sb_temp);
            failed_temp = M-num_satisfied_temp;

           if failed_init == 0
                break;
           end 

           if failed_temp < failed_init
                %fprintf("Updating failed_init: %d -> %d\n", failed_init, failed_temp);
                d_gf16 = temp_gf16;             
                failed_init = failed_temp;
            
                % Save best codeword so far
                if failed_temp < min_failed
                    min_failed = failed_temp;
                    best_d_gf16 = temp_gf16;
                    %fprintf("\nbest codeword updated!\n");    
                end
            end  
    end

    seqgf = best_d_gf16;

    % start of postprocessing     

    % Compute syndrome for the decoded sequence
    S = decod_prod(seqgf, h, CN_lst, mul_mat, add_mat);
    Sb = double(S == 0);  % 1 = satisfied, 0 = unsatisfied

    % Identify unsatisfied check nodes
    unsatisfied_indices = find(Sb == 0);

    % Get variable nodes connected to these unsatisfied check nodes
    VN_candidates = [];
    for i = 1:length(unsatisfied_indices)
        idx = unsatisfied_indices(i);
        VN_list = CN_lst{idx};
        VN_candidates = [VN_candidates, VN_list];
    end
    VN_candidates = unique(VN_candidates);  % remove duplicates

    % Start with best codeword as current decoded codeword
    best_seq = seqgf;
    best_failed = M - sum(Sb);

    for i = 1:length(VN_candidates)
        vn_idx = VN_candidates(i);
        current_symbol = seqgf(vn_idx);
        
        for test_sym = 0:15
            if test_sym ~= current_symbol
                temp_seq = seqgf;
                temp_seq(vn_idx) = test_sym;
                S_check = decod_prod(temp_seq, h, CN_lst, mul_mat, add_mat);
                Sb_check = double(S_check == 0);
                failed_now = M - sum(Sb_check);
    
                if failed_now < best_failed
                    best_failed = failed_now;
                    best_seq = temp_seq;
                end
            end
        end
    end

seqgf = best_seq;
end


