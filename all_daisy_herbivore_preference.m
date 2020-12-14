clear all
L_list = linspace(1.05, 1.5, 100)
S0 = 1380
Aw = 0.75
Ab = 0.25
Ag = 0.5
K = 0.6
b = 3.265e-3
T0 = 295.7
D = 0.3
sigma = 5.67e-8
tot_time=10000;
dt = 0.1
B1 = 0.05; % Predation rate
B2 = 0.1
F = 0.2; % Death rate for herbivore
G = 5 % Growth rate for herbivore due to grazing

Niter=ceil(tot_time/dt);
for j=1:length(L_list)
    L = L_list(j);
    Cw(1) = 0.01;
    Cb(1) = 0.01;
    x(1) = 1;
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
        Cb(i+1)=dt*(bb(i)*(1-Cb(i)-Cw(i))*Cb(i)-B1*x(i)*Cb(i))+Cb(i);
        Cw(i+1)=dt*(bw(i)*(1-Cb(i)-Cw(i))*Cw(i)-B2*x(i)*Cw(i))+Cw(i);
        x(i+1)=dt*(-F*x(i)+G*x(i)*(Cb(i+1)+Cw(i+1)))+x(i);
    end
    Cwlast(j) = Cw(end);
    Cblast(j) = Cb(end);
    Cglast(j) = Cg(end);
    xlast(j) = x(end);
end
yyaxis left
plot(L_list, Cwlast, '-rd')
hold on
plot(L_list, Cblast, '-bo')
hold on
plot(L_list, Cglast, '-gx')
title('Equilibrium populations vs. Luminosity', 'Fontsize', 30)
xlabel('Luminosity', 'Fontsize', 30)
ylabel('Daisy Surface Fractions', 'Fontsize', 30)
yyaxis right
plot(L_list, xlast, '-p')
legend('White daisies', 'Black daisies', 'Bare Ground', 'Number of herbivores', 'Fontsize', 20)
ax=gca
ax.FontSize=16
ylabel('Number of herbivores', 'Fontsize',  20)