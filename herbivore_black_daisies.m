clear all
L = 1.25;
S0 = 1380;
Cb(1) = 0.01;
Ab = 0.25;
Ag = 0.5;
K = 0.6;
b = 3.265e-3;
T0 = 295.7;
sigma = 5.67e-8
tot_time=1000;
dt = 0.1
Niter=ceil(tot_time/dt);

% Initial conditions for Lotka-Volterra model:
x(1) = 1; % 2 prey individuals per unit area
B = .05; % Predation rate
F = 0.2; % Death rate for herbivore
G = 5 % Growth rate for herbivore due to grazing

for i=1:Niter
    time(i+1) = dt*i;
    Cg(i) = 1-Cb(i);
    A(i) = Ag*Cg(i)+Ab*Cb(i);
    Ta4(i) = (L*S0*(1-A(i)))/(4*sigma);
    Ts4(i) = 2*Ta4(i);
    Tb4(i) = (1-K)*((L*S0)/(4*sigma))*(A(i)-Ab)+Ts4(i);
    Tb(i) = Tb4(i)^(1/4);
    bb(i) = 1-b*(T0-Tb(i))^2;
    if bb(i)<0
        bb(i)=0;
    end
    if x(i)<0
        x(i)=0;
    end
    if Cb(i)<0
        Cb(i)=0;
    end
    Cb(i+1)=dt*(bb(i)*(1-Cb(i))*Cb(i)-B*x(i)*Cb(i))+Cb(i);
    x(i+1)=dt*(-F*x(i)+G*x(i)*(Cb(i+1)))+x(i);
    
end
% ta = (ta4)^(1/4)
subplot(2, 2, 1)
yyaxis left
plot(time, Cb, 'k')
% title('Fraction covered by black and white daisies')
xlabel('Time')
ylabel('Fraction')
hold off
% legend(' Fraction covered by White daisies', 'Fraction covered by Black daisies')
yyaxis right
plot(time, x, 'c')
% title('Populations of predator and prey per unit area')
xlabel('Time')
ylabel('Number per unit area')
legend('Fraction covered by Black daisies', 'Herbivores per unit area')
%xlim([0 20])

subplot(2, 2, 2)
plot([time(1:end-1)], Tb, 'k')
title('Temperature')
xlabel('Time')
ylabel('Temp')
legend('Black daisies Temp')

subplot(2, 2, 4)
plot(x, Cb, 'k')
xlabel('Herbivores per unit area')
ylabel('Fraction covered by black daisies')
hold on
% quiver(x,Cb,gradient(x),gradient(Cb))
title('Herbivore/black daisy phase space')