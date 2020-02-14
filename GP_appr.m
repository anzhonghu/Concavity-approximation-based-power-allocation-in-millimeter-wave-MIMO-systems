%GP_appr
gammat = 1.43 *ones(K, 1);
alphat = 1e3;
betat = 1e-1;
pkt = zeros(K, 1);
pkt1 = real(rho0) * ones(K, 1);
t_count = 0;
while norm(pkt1-pkt) > 1e-2 * mean(pkt)
    t_count = t_count + 1;
    if t_count > 20
        break;
    else
    end
    deltapk1 = norm(pkt1 - pkt);
    pkt = pkt1;
    lambdat = K / (log(2) * rho);
    flag_lambda = 1;
    pos_c = 0;
    neg_c = 0;
    while abs(flag_lambda-1) < 1e-2
        for k = 1 : K
            pksum = 0;
            for kk = 1 : K
                if kk == k
                else
                    pksum = pksum + abs(akkps(kk, k)) * gammat(kk, 1);
                end
            end
            pkt1(k, 1) =  1 / log(2) / (lambdat * bkks(k, 1) + pksum);
        end
        pksum = 0;
        for kk = 1 : K
            pksum = pksum + bkks(kk, 1) * pkt1(kk, 1);
        end
        if pksum - rho > 1e-2 * rho
            lambdat = lambdat * 1.1;
            pos_c = pos_c + 1;
        else
            if pksum - rho < -1e-2 * rho
                lambdat = lambdat * 0.91;
                neg_c = neg_c + 1;
            else
                flag_lambda = 0;
            end
        end
        if pos_c > 10 || neg_c > 10
            break;
        else
        end
    end
    for k = 1 : K
        ezkt = 1 / log(2) / gammat(k, 1) - 1;
        pksum = 0;
        for kk = 1 : K
            if kk == k
            else
                pksum = pksum + abs(akkps(k, kk)) * pkt(kk, 1);
            end
        end
        gammat(k, 1) = max(gammat(k, 1) + real(betat * (-pksum + ezkt)), 1e-3);
    end
end
pkt = pkt1;
pksum = 0;
for kk = 1 : K
    pksum = pksum + bkks(k, 1) * pkt(kk, 1);
end
pkt = pkt * rho / pksum;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1 : K
    interf = 0;
    for kk = 1 : K
        if kk == k
        else
            interf = interf + abs(akkps(k, kk)) * pkt(kk, 1);
        end
    end
    interf = interf + 1;
    signalpow = abs(akks(k, 1)) * pkt(k, 1);
    sinr = real(signalpow  / interf);
    capacity(signalpo_n, 1) = capacity(signalpo_n, 1) + log2(1 + sinr);
end
