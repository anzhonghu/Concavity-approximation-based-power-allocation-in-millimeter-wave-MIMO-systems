
D3 = zeros(Mb, 1);
diff = 1e4 * ones(n_arr(3, 1), 1);
H3_e =  U3' * (H3);
[~, uhe3in] = sort(abs(H3_e), 'descend');
for k = 1 : K
    counttemp = 0;
    for n = 1 : n_arr(3, 1)
        if diff(uhe3in(n,k), 1) > 0
            diff(uhe3in(n,k), 1) = 0;
            counttemp = counttemp + 1;
            D3((k-1)*beam_numberforMS+counttemp, 1) = uhe3in(n,k);
        else
        end
        if counttemp < beam_numberforMS
        else
            break;
        end
    end
end
H3eb = H3_e(D3, :);%%Mb*K
U3b = U3(:, D3);
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%MMSE uspa%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%Initialization%%%%%%%%%%%%%%%%
%%%%%%%%%Initialization of algorithm 3%%%%%%%%%%%%%
[Vh3eb, H3d, ~] = svd(H3eb' * H3eb);
Kq = K;
AA = Vh3eb * (eye(Kq*Mm) - inv(eye(Kq*Mm) + rho/Kq/Mm * H3d)) * Vh3eb';
BB = Vh3eb * H3d /(Kq*Mm/rho*eye(Kq*Mm) + H3d)^2 * Vh3eb';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%algorithm 2%%%%%%%%%%%%%%%%%
indH3d = ones(Kq, 1);
D3q = D3;
H3q = H3_e;
akks = zeros(Kq, 1);
bkks = zeros(Kq, 1);
for k = 1 : Kq
    for m1=1:Mm
        for m2=1:Mm
            akks(k, 1) = akks(k, 1) + AA((k-1)*Mm+m1, (k-1)*Mm+m2);
        end
    end
    for m1=1:Mm
        for m2=1:Mm
            bkks(k, 1) = bkks(k, 1) + BB((k-1)*Mm+m1, (k-1)*Mm+m2);
        end
    end
end
akks = real(akks);
bkks = real(bkks);
akkps = zeros(Kq, Kq);
for k = 1 : Kq
    for kk = 1 : Kq
        if k == kk
            akkps(k, kk) = akks(k, 1);
        else
            for m1=1:Mm
                for m2=1:Mm
                    akkps(k, kk) = akkps(k, kk) + AA((k-1)*Mm+m1, (kk-1)*Mm+m2);
                end
            end
        end
    end
end
akkps = real(akkps);
akks = akks.^2;
akkps = akkps.^2;
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%MMSE%%%%%%%%%%%%%%%%%%%%
rho0 = rho / sum(bkks);
pks = rho0 * ones(Kq, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
partid = zeros(Kq, 1);
partiddb = zeros(Kq, 1);
patimean = 0;
akksmean = mean(akks(k, 1));
sumbkks = 0;

for k = 1 : Kq
    gk =  1;
    for kk = 1 : Kq
        if k == kk
        else
            gk = gk + abs(akkps(k, kk)) * pks(kk, 1);
        end
    end
    gk = akks(k, 1) * pks(k, 1) / gk;
    if k < Kq
        k0 = k + 1;
    else
        k0 = k - 1;
    end
    xk = 0;
    for kk = 1 : Kq
        if kk == k0
        else
            xk = xk + bkks(kk, 1) / bkks(k0, 1) * pks(kk, 1);%b~kk * pkk
        end
    end
    xk = rho / bkks(k0, 1)  - xk;
    xk = abs(akkps(k, k0)) * xk + 1;%akk0*()+1
    for kk = 1 : Kq
        if kk == k0 || kk == k
        else
            xk = xk + abs(akkps(k, kk)) * pks(kk, 1);%akkk*pkkk
        end
    end
    akk0 = abs(akkps(k, k0));
    bqk = bkks(k, 1) / bkks(k0, 1);
    partid(k, 1) = 1 / (1+gk) * akks(k, 1) * (xk + akk0 * bqk * pks(k, 1)) / (xk^2);
    patimean = patimean + partid(k, 1) / bkks(k, 1);
    partiddb(k, 1) = partid(k, 1) / bkks(k, 1);
end
patimean = patimean / Kq;
patidevi = 0;
patidevi_s = 1e10;
for k = 1 : Kq
    patidevi = patidevi +  (partiddb(k, 1) - patimean)^2;
end
patidevi_comp = min(partiddb);
ite_c = 0;
partiddb_s = 0;
flag_pp = zeros(K, 1);
% ky_flag = zeros(K,1);

patidevi_st(1,k_ind) = patidevi;
while patidevi > patidevi_comp * 1e-12
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    partiddb_sss = partiddb;
    patidevi_sss = patidevi;
    ite_c = ite_c + 1;
    if ite_c > 10*K
        break;
    else
    end
    [~, kx] = max(partiddb);
    for k = 1 : K
        if indH3d(k, 1) < 1e-5
            partiddb(k, 1) = 1e10;
        else
        end
    end
    [~, ky] = min(partiddb);
    %%%%%%%%%%%%%%%%%
    gkx =  1;
    for kk = 1 : K
        if kk == kx
        else
            gkx = gkx + abs(akkps(kx, kk)) * pks(kk, 1);
        end
    end
    gkx = akks(kx, 1) * pks(kx, 1) / gkx;
    k0 = 1;
    while k0==kx || k0==ky
        k0 = k0 + 1;
    end
    if k0 > K
        break;
    else
    end
    xkx = 0;
    for kk = 1 : K
        if kk == k0
        else
            xkx = xkx + bkks(kk, 1) / bkks(k0, 1) * pks(kk, 1);%b~kk * pkk
        end
    end
    xkx = rho / bkks(k0, 1)  - xkx;
    xkx = abs(akkps(kx, k0)) * xkx + 1;%akk0*()+1
    for kk = 1 : K
        if kk == k0 || kk == kx
        else
            xkx = xkx + abs(akkps(kx, kk)) * pks(kk, 1);%akkk*pkkk
        end
    end
    akxk0 = abs(akkps(kx, k0));
    bqkx = bkks(kx, 1) / bkks(k0, 1);
    xqkx = (akks(kx, 1) * (xkx + akxk0 * bqkx * pks(kx, 1)));
    parti2pkx = -1/ (1+gkx)^2 * xqkx^2 / xkx^4 + 2 / (1+gkx) * akxk0 * bqkx * (akks(kx, 1) * (xkx + akxk0 * bqkx * pks(kx, 1))) / xkx^3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gky =  1;
    for kk = 1 : K
        if kk == ky
        else
            gky = gky + abs(akkps(ky, kk))* pks(kk, 1);
        end
    end
    gky = akks(ky, 1) * pks(ky, 1) / gky;
    xky = 0;
    for kk = 1 : K
        if kk == k0
        else
            xky = xky + bkks(kk, 1) / bkks(k0, 1) * pks(kk, 1);%b~kk * pkk
        end
    end
    xky = rho / bkks(k0, 1)  - xky;
    xky = abs(akkps(ky, k0)) * xky + 1;%akk0*()+1
    for kk = 1 : K
        if kk == k0 || kk == ky
        else
            xky = xky + abs(akkps(ky, kk)) * pks(kk, 1);%akkk*pkkk
        end
    end
    akyk0 = abs(akkps(ky, k0));
    bqky = bkks(ky, 1) / bkks(k0, 1);
    xqky = (akks(ky, 1) * (xky + akyk0 * bqky * pks(ky, 1)));
    parti2pky = -1/ (1+gky)^2 * xqky^2 / xky^4 + 2 / (1+gky) * akyk0 * bqky * (akks(ky, 1) * (xky + akyk0 * bqky * pks(ky, 1))) / xky^3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    akxky = abs(akkps(kx, ky));
    xqkx_modi = akks(kx, 1)*((1+gkx)*(xkx+2*akxk0*bqkx * pks(kx, 1))-xqkx/xkx*pks(kx, 1));
    parti2fkx_pkxpky =  xqkx_modi * (akxk0 * bqky - akxky) / (1+gkx)^2  / xkx^3;
    %%%%%%%%%%%%%%%%%%%%%%
    akykx = abs(akkps(ky, kx));
    xqky_modi = akks(ky, 1)*((1+gky)*(xky+2*akyk0*bqky * pks(ky, 1))-xqky/xky*pks(ky, 1));
    parti2fky_pkxpky =  xqky_modi * (akyk0 * bqkx - akykx) / (1+gky)^2  / xky^3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bky = bkks(ky, 1);
    bkx = bkks(kx, 1);
    deltapkx = (partiddb(kx, 1)-partiddb(ky, 1)) / (parti2fkx_pkxpky/bky + parti2fkx_pkxpky/bky - parti2pky * bkx / bky^2 - parti2pkx / bkx);
    deltapky = -deltapkx * bkx / bky;
    %%%%%%%%%%%%%%%%%%%%%%%%%
    pksx_sss = pks(kx, 1);
    pksy_sss = pks(ky, 1);
    if pks(ky, 1) + deltapky <= 0
        pks(kx, 1) = pks(kx, 1) + pks(ky, 1) * bky / bkx;
        pks(ky, 1) = 0;
        indH3d(ky, 1) = 0;
    else
        pks(kx, 1) = pks(kx, 1) + deltapkx;
        pks(ky, 1) = pks(ky, 1) + deltapky;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Kqq = sum(indH3d);
    if Kqq < 3
        break;
    else
    end
    Kq = Kqq;
    partid = zeros(K, 1);
    partiddb = zeros(K, 1);
    patimean = 0;
    for k = 1 : K
        if indH3d(k, 1) < 1e-5
        else
            gk =  1;
            for kk = 1 : K
                if kk == k
                else
                    gk = gk + abs(akkps(k, kk)) * pks(kk, 1);
                end
            end
            gk = akks(k, 1) * pks(k, 1) / gk;
            if k < K
                k0 = k + 1;
            else
                k0 = k - 1;
            end
            xk = 0;
            for kk = 1 : K
                if kk == k0
                else
                    xk = xk + bkks(kk, 1) / bkks(k0, 1) * pks(kk, 1);%b~kk * pkk
                end
            end
            xk = rho / bkks(k0, 1)  - xk;
            xk = abs(akkps(k, k0)) * xk + 1;%akk0*()+1
            for kk = 1 : Kq
                if kk == k0 || kk == k
                else
                    xk = xk + abs(akkps(k, kk)) * pks(kk, 1);%akkk*pkkk
                end
            end
            akk0 = abs(akkps(k, k0));
            bqk = bkks(k, 1) / bkks(k0, 1);
            partid(k, 1) = 1 / (1+gk) * akks(k, 1) * (xk + akk0 * bqk * pks(k, 1)) / (xk^2);
            patimean = patimean + partid(k, 1) / bkks(k, 1);
            partiddb(k, 1) = partid(k, 1) / bkks(k, 1);
        end
    end
    partiddb_s = partiddb(ky, 1);
    patimean = patimean / Kq;
    patidevi = 0;
    for k = 1 : K
        if indH3d(k, 1) < 1e-5
        else
            patidevi = patidevi +  (partiddb(k, 1) - patimean)^2;
        end
    end
    if patidevi_s *10 < patidevi
        pks(kx, 1) = pksx_sss;
        pks(ky, 1) = pksy_sss;
        partiddb = partiddb_sss;
        patidevi = patidevi_sss;
        ite_c = ite_c - 1;
        break;
    else
    end
    if patidevi_s <= patidevi && 0 ~= indH3d(ky, 1)%%&& 1 ~= flag_pp(ky,1)
        break;
    else
    end
    patidevi_s = patidevi;
    partiddb_comp = partiddb;
    for k = 1 : K
        if indH3d(k, 1) < 1e-5
            partiddb_comp(k, 1) = 1e10;
        else
        end
    end
    [patidevi_comp, ~] = min(partiddb_comp);
    patidevi_st(ite_c+1,k_ind) = patidevi;
end
if 1==k_ind
ite_c_s = ite_c;
else
    if ite_c>ite_c_s
        ite_c_s = ite_c;
    else
    end
end