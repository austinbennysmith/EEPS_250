clear all
L = 1.55;
S0 = 1380;
Cw(1) = 0.01;
Cb(1) = 0.01;
Aw = 0.75;
Ab = 0.25;
Ag = 0.5;
K = 0.6;
b = 3.265e-3;
T0 = 295.7;
sigma = 5.67e-8
tot_time=500;
dt = 0.1
Niter=ceil(tot_time/dt);

% Initial conditions for Lotka-Volterra model:
x(1) = 1; % number of individual herbivores
B = 0.05; % Predation rate
F = 0.2; % Death rate for herbivore
G = 5 % Growth rate for herbivore due to grazing

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
     if x(i)<0
         x(i)=0;
     end
     if Cb(i)<0
         Cb(i)=0;
     end
     if Cw(i)<0
         Cw(i)=0;
     end
    Cb(i+1)=dt*(bb(i)*(1-Cb(i)-Cw(i))*Cb(i)-B*x(i)*Cb(i))+Cb(i);
    Cw(i+1)=dt*(bw(i)*(1-Cb(i)-Cw(i))*Cw(i)-B*x(i)*Cw(i))+Cw(i);
    x(i+1)=dt*(-F*x(i)+G*x(i)*(Cb(i+1)+Cw(i+1)))+x(i);
    
end
% ta = (ta4)^(1/4)
subplot(2, 2, 1)
yyaxis left
plot(time, Cw, 'm')
hold on
plot(time, Cb, 'k')
% title('Fraction covered by black and white daisies')
xlabel('Time', 'Fontsize', 30)
ylabel('Fraction', 'Fontsize', 30)
hold off
% legend(' Fraction covered by White daisies', 'Fraction covered by Black daisies')
yyaxis right
plot(time, x, 'c')
% title('Populations of predator and prey per unit area')
xlabel('Time',  'Fontsize', 30)
ylabel('Number of herbivores', 'Fontsize', 30)
legend(' Fraction covered by White daisies', 'Fraction covered by Black daisies', 'Number of Herbivores', 'Fontsize', 15)
title('Daisy fraction/number of herbivores vs time', 'Fontsize', 30)
%xlim([0 20])

subplot(2, 2, 2)
plot([time(1:end-1)], Tw, 'm')
hold on
plot([time(1:end-1)], Tb, 'k')
title('Temperature', 'Fontsize', 30)
xlabel('Time', 'Fontsize', 30)
ylabel('Temp (K)', 'Fontsize', 30)
legend('White daisies Temp', 'Black daisies Temp', 'Fontsize', 15)

subplot(2, 2, 3)
plot(x, Cw, 'm')
xlabel('Number of Herbivores', 'Fontsize', 30)
ylabel('Fraction covered by wite daisies', 'Fontsize', 30)
title('Herbivore/white daisy phase space', 'Fontsize', 30)
hold on
% quiver(x,Cw,gradient(x),gradient(Cw))

subplot(2, 2, 4)
plot(x, Cb, 'k')
xlabel('Number of Herbivores', 'Fontsize', 30)
ylabel('Fraction covered by black daisies', 'Fontsize', 30)
hold on
% quiver(x,Cb,gradient(x),gradient(Cb))
title('Herbivore/black daisy phase space', 'Fontsize',  30)