%SCA_appr
pkt = rho0 * ones(K, 1);
ht = 1;
ht1 = 0;
alphaback = 0.5;
betaback = 0.5;
t_count = 0;
while abs(ht1-ht) > 1e-2 * abs(ht)
    t_count = t_count + 1;
    if t_count > 20
        break;
    else
    end
    ht = ht1;
    lambdasca = 1e3;
    patif = zeros(K, 1);
    for k0 = 1 : K
        for k =1 : K
            nominat = 0;
            for kk = 1 : K
                nominat = nominat + abs(akkps(k, kk)) * pkt(kk, 1);
            end
            patif(k0, 1) = patif(k0, 1) + 1 / log(2) / (nominat+1) * abs(akkps(k, k0));
        end
    end
    patif = - patif;
    patig = zeros(K, 1);
    for k0 = 1 : K
        for k =1 : K
            nominat = 0;
            for kk = 1 : K
                if kk == k
                else
                    nominat = nominat + abs(akkps(k, kk)) * pkt(kk, 1);
                end
            end
            if k0 == k
            else
                patig(k0, 1) = patig(k0, 1) + 1 / log(2) / (nominat+1) * abs(akkps(k, k0));
            end
        end
    end
    patig = - patig;
    patih = patif - patig;
    patih = patih - mean(patih);
    deltap = -patih ./ bkks * mean(bkks);
    %%%%%%%%%%%%%%%%%%%%%%%%
    t = 1;
    fxp = 1e10;
    fpt = 0;
    gpt = 0;
    for k =1 : K
        nominat = 0;
        for kk = 1 : K
            nominat = nominat + abs(akkps(k, kk)) * pkt(kk, 1);
        end
        fpt = fpt + log2(nominat+1);
    end
    fpt = -fpt;
    for k =1 : K
        nominat = 0;
        for kk = 1 : K
            if kk == k
            else
                nominat = nominat + abs(akkps(k, kk)) * pkt(kk, 1);
            end
        end
        gpt = gpt + log2(nominat+1);
    end
    gpt = -gpt;
    fx = fpt - gpt;
    while fxp - fx - alphaback * t * patih' * deltap > 0
        t = betaback * t;
        pktemp = pkt + deltap * t;
        if min(pktemp) < 0
            continue;
        else
        end
        if t < 1e-6
            break;%%%function h is too flat, and is actually the case in simulation, which means average power allocation
        else
        end
        fpt = 0;
        for k =1 : K
            nominat = 0;
            for kk = 1 : K
                nominat = nominat + abs(akkps(k, kk)) * pktemp(kk, 1);
            end
            fpt = fpt + log2(nominat+1);
        end
        fpt = -fpt;
        fxp = fpt - gpt - patig' * deltap * t;
    end
    pbar = pkt + deltap * t;
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    fpt = 0;
    for k =1 : K
        nominat = 0;
        for kk = 1 : K
            nominat = nominat + abs(akkps(k, kk)) * pbar(kk, 1);
        end
        fpt = fpt + log2(nominat+1);
    end
    fpt = -fpt;
    ht1 = fpt - gpt - patig' * (pbar - pkt);
    %%%%%%%%%%%%%%%%%%%5
    pkt = pbar;
    if t < 1e-6
        break;
    else
    end
end
pkt = real(pkt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1 : K
    akkp = 0;
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
    capacity(signalpo_n, 2) = capacity(signalpo_n, 2) + log2(1 + sinr);
end