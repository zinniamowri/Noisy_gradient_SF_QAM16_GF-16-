function [G, K, pivotCols, freeCols, H_rref] = gf16_generator_from_H(H, add_mat, mul_mat, div_mat)
%GF16_GENERATOR_FROM_H Build generator matrix G from parity-check matrix H over GF(16).
%
% Inputs:
%   H       : MxN, entries in {0..15} (0 means no edge / zero element)
%   add_mat : 16x16 lookup table for GF add, indexed by (a+1,b+1)
%   mul_mat : 16x16 lookup table for GF multiply
%   div_mat : 16x16 lookup table for GF divide: a/b (b!=0)
%
% Outputs:
%   G        : KxN generator matrix over GF(16) (entries 0..15), satisfies H*G' = 0
%   K        : dimension = N - rank(H)
%   pivotCols: pivot column indices in RREF
%   freeCols : free column indices (these correspond to message symbol positions)
%   H_rref   : RREF form of H (over GF(16)) used for construction, row-reduced echelon form

    H = uint8(H);
    [M, N] = size(H);

   
    H_rref = H;
    pivotCols = zeros(1, min(M,N));
    r = 0;  % rank
    row = 1;

    for col = 1:N
        if row > M
            break;
        end

        % Find pivot row with nonzero in this column
        piv = 0;
        for rr = row:M
            if H_rref(rr,col) ~= 0
                piv = rr;
                break;
            end
        end
        if piv == 0
            continue; 
        end

       
        if piv ~= row
            tmp = H_rref(row,:);
            H_rref(row,:) = H_rref(piv,:);
            H_rref(piv,:) = tmp;
        end

        
        pivotVal = H_rref(row,col); % nonzero
        invPivot = gf_div(uint8(1), pivotVal, div_mat); 
        H_rref(row,:) = gf_row_scale(H_rref(row,:), invPivot, mul_mat);

        % Eliminate this column in ALL other rows
        for rr = 1:M
            if rr == row
                continue;
            end
            factor = H_rref(rr,col);
            if factor ~= 0
                
                H_rref(rr,:) = gf_row_add_scaled(H_rref(rr,:), H_rref(row,:), factor, add_mat, mul_mat);
            end
        end

        r = r + 1;
        pivotCols(r) = col;
        row = row + 1;
    end

    pivotCols = pivotCols(1:r);
    allCols = 1:N;
    freeCols = setdiff(allCols, pivotCols, 'stable');
    K = numel(freeCols);


    G = zeros(K, N, 'uint8');

    for t = 1:K
        x = zeros(1, N, 'uint8');
        fcol = freeCols(t);
        x(fcol) = 1;

        
        for i = 1:r
            pcol = pivotCols(i);

            acc = uint8(0);
            
            for jj = 1:K
                jcol = freeCols(jj);
                hij = H_rref(i, jcol);
                if hij ~= 0 && x(jcol) ~= 0
                    acc = gf_add(acc, gf_mul(hij, x(jcol), mul_mat), add_mat);
                end
            end

            
            x(pcol) = acc;
        end

        G(t,:) = x;
    end
end


function c = gf_add(a,b,add_mat)
    c = uint8(add_mat(double(a)+1, double(b)+1));
end

function c = gf_mul(a,b,mul_mat)
    c = uint8(mul_mat(double(a)+1, double(b)+1));
end

function c = gf_div(a,b,div_mat)
    
    c = uint8(div_mat(double(a)+1, double(b)+1));
end

function row2 = gf_row_scale(row, scalar, mul_mat)
    row2 = row;
    if scalar == 0
        row2(:) = 0;
        return;
    end
    for k = 1:numel(row)
        if row(k) ~= 0
            row2(k) = gf_mul(row(k), scalar, mul_mat);
        else
            row2(k) = 0;
        end
    end
end

function rowOut = gf_row_add_scaled(rowA, rowPivot, factor, add_mat, mul_mat)
    
    rowOut = rowA;
    for k = 1:numel(rowA)
        term = uint8(0);
        if rowPivot(k) ~= 0 && factor ~= 0
            term = gf_mul(factor, rowPivot(k), mul_mat);
        end
        rowOut(k) = gf_add(rowA(k), term, add_mat);
    end
end
