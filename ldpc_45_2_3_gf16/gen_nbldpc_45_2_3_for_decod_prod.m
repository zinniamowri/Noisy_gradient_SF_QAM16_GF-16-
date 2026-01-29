function [h, str_cn_vn, Hb] = gen_nbldpc_45_2_3_for_decod_prod(seed, maxTries)

    if nargin < 1 || isempty(seed), seed = 1; end
    if nargin < 2 || isempty(maxTries), maxTries = 20000; end
    rng(seed);

    N = 45; dv = 2; dc = 3;
    M = (N*dv)/dc;

    Hb = false(M,N);

    for attempt = 1:maxTries
        Hb(:) = false;
        cn_deg = zeros(M,1);
        commonCount = zeros(N,N,'uint8');
        ok = true;

        for v = randperm(N)
            chosen = [];
            for e = 1:dv
                cand = find(cn_deg < dc);
                cand = setdiff(cand, chosen);
                if isempty(cand), ok=false; break; end

                mask = false(size(cand));
                for kk = 1:numel(cand)
                    c = cand(kk);
                    nbr = find(Hb(c,:));
                    if isempty(nbr) || all(commonCount(nbr,v) == 0)
                        mask(kk) = true;
                    end
                end

                validCand = cand(mask);    % <<< ADD THIS LINE
                
                if isempty(validCand)
                    ok = false;
                    break;
                end
                
                d = cn_deg(validCand);
                best = validCand(d == min(d));
                cPick = best(randi(numel(best)));


                Hb(cPick,v) = true;
                cn_deg(cPick) = cn_deg(cPick) + 1;

                nbr = find(Hb(cPick,:));
                nbr(nbr==v) = [];
                for u = nbr
                    commonCount(u,v)=commonCount(u,v)+1;
                    commonCount(v,u)=commonCount(u,v);
                    if commonCount(u,v)>1
                        ok=false; break;
                    end
                end
                if ~ok, break; end
                chosen(end+1)=cPick;
            end
            if ~ok, break; end
        end

        if ok && all(cn_deg==dc)
            break;
        end
    end

    assert(ok,'Failed to generate H without 4-cycles');

    % Assign GF(16) weights
    h = zeros(M,N,'uint8');
    [r,c] = find(Hb);
    w = randi([1 15],numel(r),1,'uint8');
    for k=1:numel(r)
        h(r(k),c(k)) = w(k);
    end

    % Build CN list (EXACTLY what decod_prod expects)
    str_cn_vn = cell(M,1);
    for i=1:M
        str_cn_vn{i} = find(h(i,:)~=0);
    end
end
