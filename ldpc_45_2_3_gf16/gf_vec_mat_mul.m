function c = gf_vec_mat_mul(u, G, add_mat, mul_mat)
% u : 1xK message over GF(16)
% G : KxN generator matrix over GF(16)
% c : 1xN codeword over GF(16)

[K, N] = size(G);
c = zeros(1, N, 'uint8');

for j = 1:N
    acc = uint8(0);
    for k = 1:K
        if u(k) ~= 0 && G(k,j) ~= 0
            acc = uint8(add_mat(double(acc)+1, ...
                double(mul_mat(double(u(k))+1, double(G(k,j))+1))+1));
        end
    end
    c(j) = acc;
end
end
