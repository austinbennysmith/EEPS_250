clear all
L_list = linspace(0.7, 2, 100)

S0 = 1380
Aw = 0.75
Ab = 0.25
Ag = 0.5
K = 0.6
b = 3.265e-3
T0 = 295.7
D = 0.3
sigma = 5.67e-8
tot_time=1000;
dt = 0.1

Niter=ceil(tot_time/dt);
for j=1:length(L_list)
    L = L_list(j);
    Cw(1) = 0.01;
    Cb(1) = 0.01;
    for i=1:Niter
        time(i+1) = dt*i;
        Cg(i) = 1-Cw(i)-Cb(i);
        A(i) = Aw*Cw(i)+Ag*Cg(i)+Ab*Cb(i);
        Ta4(i) = (L*S0*(1-A(i)))/(4*sigma);
        Ts4(i) = 2*Ta4(i);
        Tw4(i) = (1-K)*((L*S0)/(4*sigma))*(A(i)-Aw)+Ts4(i);
        Tw(i) = Tw4(i)^(1/4);
        Tb4(i) = (1-K)*((L*S0)/(4*sigma))*(A(i)-Ab)+Ts4(i);
        Tb(i) = Tb4(i)^(1/4);
        bb(i) = 1-b*(T0-Tb(i))^2;
        bw(i) = 1-b*(T0-Tw(i))^2;
        if bw(i) < 0
            bw(i)=0;
        end
        if bb(i)<0
            bb(i)=0;
        end
        Cb(i+1) = Cb(i)*(1+dt*(bb(i)*(1-Cb(i)-Cw(i))-D));
        Cw(i+1) = Cw(i)*(1+dt*(bw(i)*(1-Cb(i)-Cw(i))-D));
    end
    Cwlast(j) = Cw(end);
    Cblast(j) = Cb(end);
    Cglast(j) = Cg(end);
    bblast(j) = bb(end)
    bwlast(j) = bw(end)
    Tslast(j) = Ts4(end)^(1/4)
end
% subplot(1, 2, 1)
plot(L_list, Cwlast, '-rd')
hold on
plot(L_list, Cblast, '-bo')
hold on
plot(L_list, Cglast, '-gx')
title('Equilibrium Surface Fractions vs. Luminosity', 'FontSize', 40)
xlabel('Luminosity', 'Fontsize', 40)
ylabel('Surface Fractions', 'FontSize', 40)
legend('White daisies', 'Black daisies', 'Bare Ground')
% subplot(1, 2, 2)
% yyaxis left
% plot(L_list, bwlast, '-r', 'Linewidth', 2)
% hold on
% plot(L_list, bblast, '-b', 'Linewidth', 2)
% xlabel('Luminosity')
% ylabel('Birth rate')
% title('Equilibrium state birth rates as functions of luminosity')
% yyaxis right
% plot(L_list, Tslast, '-k', 'Linewidth', 2)
% legend('White daisies', 'black daisies', 'Equillibrium surface temp')