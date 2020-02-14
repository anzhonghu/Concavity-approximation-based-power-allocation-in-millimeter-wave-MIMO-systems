clear;
close all;
naz = [7; 16; 32];
nel = [36; 76; 151];
n_arr = naz .* nel;

Rmin = 10;
Rmax = 100;
rho_store = [1e6;10^(6.5);1e7;10^(7.5);1e8;];%[1e9; 10^(9.5); 1e10; 10^(10.5); 1e11;]*1e-2;
height = 10;
f = 80*10^9;%1G bandwidth
lambda = 3 * 10^8 / f;
miu = 0.5;
Spad = 5;
Mm = 4;
Pp = 2;
pblock = 0.1;
betam_min = atan(height/Rmax);
betam_max = atan(height/Rmin);
% beta_m = 0.5*(betam_min+betam_max);%larger than 0.5*(betam_min+betam_max)
beta_m = 10.3 * pi / 180;%-phi_m in the EL paper
Nite = 1e2;


capacity = zeros(length(rho_store),  4);
angle = zeros(1, 2);
U3 = zeros(n_arr(3, 1), n_arr(3, 1));
for nx = 1 : naz(3, 1)
    for ny = 1 : nel(3, 1)
        angle(1, 2) = (-1+2*ny/nel(3, 1));%el
        angle(1, 1) = (-1+2*nx/naz(3, 1));%az
        n = (ny - 1) * naz(3, 1) + nx;
        for mx = 0 : naz(3, 1)-1
            for my = 0 : nel(3, 1)-1
                m = my * naz(3, 1) + 1 + mx;
                U3(m, n) = exp(-1i * 2 * pi * miu * ((mx-0.5*(naz(3, 1)-1)) * angle(1, 1) + (my-0.5*(nel(3, 1)-1)) * angle(1, 2))) / sqrt(n_arr(3, 1));
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_ar_in =  3;
patidevi_st = zeros(1e3, 5);
K_store = [10; 30; 50; 70; 100];
for k_ind = 1 : 5
    K = K_store(k_ind, 1);
    Mb = 2 * K;
    el_in = zeros(Mb, 1);
    az_in = zeros(Mb, 1);
    beam_numberforMS = floor(Mb/K);
    c_km = zeros(K, beam_numberforMS);
    H3 = zeros(n_arr(3, 1), K*Mm);
    H3eb = zeros(Mb, K*Mm);
    rho = rho_store(length(rho_store), 1);
    pos = zeros(K, 3);
    theta = zeros(K, 1);
    phi = zeros(K, 1);
    spatail_frequ = zeros(K, 2);
    beta = zeros(K, 2);
    beam_in = zeros(Mb, 1);
    for k = 1 : K
        pos_temp = zeros(1, 2);
        while norm(pos_temp) < Rmin || norm(pos_temp) > Rmax || abs(atan(pos_temp(1, 2) / pos_temp(1, 1))) > pi / 3
            pos_temp(1, 1) = rand(1, 1) * Rmax;
            pos_temp(1, 2) = (rand(1, 1) * 2 - 1) * Rmax;
        end
        pos(k, 1:2) = pos_temp;
        pos(k, 3) = norm(pos_temp);
        pos(k, 3) = norm([pos(k, 3), height]);%distance
        phi(k, 1) = asin(pos_temp(1, 2) / sqrt(pos_temp(1, 2)^2 + (pos_temp(1, 1) * cos(beta_m) + height * sin(beta_m))^2));%az
        theta(k, 1) = asin((pos_temp(1, 1) * sin(beta_m) - height * cos(beta_m)) / pos(k, 3));%el
        spatail_frequ(k, 1) = cos(theta(k, 1)) * sin(phi(k, 1));
        spatail_frequ(k, 2) = sin(theta(k, 1));
    end
    for k = 1 : K
        phikp = rand(1,Pp+1) * 2 * pi;
        for p = 1 : Pp+1
            beta(k, 1) = sqrt(lambda^2 / (16 * pi^2 * pos(k, 3)^2)) * exp(1i * rand(1,1) * 2 * pi);
            h = zeros(n_arr(n_ar_in, 1), 1);
            if 1 == p && rand(1,1)<pblock
                continue;
            else
            end
            if 1 == p
                for nx = 0 : naz(n_ar_in, 1)-1
                    for ny = 0 : nel(n_ar_in, 1)-1
                        n = ny * naz(n_ar_in, 1) + 1 + nx;
                        h(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz(n_ar_in, 1)-1)) * sin(phi(k, 1)) * cos(theta(k, 1)) + (ny-0.5*(nel(n_ar_in, 1)-1)) * sin(theta(k, 1))));
                    end
                end
                for m = 1 : Mm
                    H3(:, (k-1)*Mm+m) = H3(:, (k-1)*Mm+m) + h * beta(k, 1) * exp(1i*pi*(m-1)*cos(phikp(1,p)));
                end
            else
                phit = phi(k, 1) +(rand(1,1)-0.5)*20*pi/180;
                thetat = theta(k, 1) +(rand(1,1)-0.5)*6*pi/180;
                for nx = 0 : naz(n_ar_in, 1)-1
                    for ny = 0 : nel(n_ar_in, 1)-1
                        n = ny * naz(n_ar_in, 1) + 1 + nx;
                        h(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz(n_ar_in, 1)-1)) * sin(phit) * cos(thetat) + (ny-0.5*(nel(n_ar_in, 1)-1)) * sin(thetat)));
                    end
                end
                for m = 1 : Mm
                    H3(:, (k-1)*Mm+m) = H3(:, (k-1)*Mm+m) + h * beta(k, 1) * exp(1i*pi*(m-1)*cos(phikp(1,p)))* 10^((-rand(1,1) * 5 - 15)*0.05);
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    
    
    cac_f_s;
    
end

ite_c = ite_c_s;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
set(h,'PaperType','A4');
axes('FontSize',16);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
semilogy(0:ite_c, patidevi_st(1:ite_c+1, :), 'k-','LineWidth',2,'MarkerSize',10)
ylim([1e-19, 1e-12])
xlabel('Number of iterations','Fontsize',20,'Fontname','Times')
ylabel('Variance of derivative','Fontsize',20,'Fontname','Times')
grid on
print(h,'-dpdf','convergence')