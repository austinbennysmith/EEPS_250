clear all
L = 1.6;
S0 = 1380;
Cw(1) = 0.01;
Aw = 0.75;
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
    Cg(i) = 1-Cw(i);
    A(i) = Ag*Cg(i)+Aw*Cw(i);
    Ta4(i) = (L*S0*(1-A(i)))/(4*sigma);
    Ts4(i) = 2*Ta4(i);
    Tw4(i) = (1-K)*((L*S0)/(4*sigma))*(A(i)-Aw)+Ts4(i);
    Tw(i) = Tw4(i)^(1/4);
    bw(i) = 1-b*(T0-Tw(i))^2;
    if bw(i)<0
        bw(i)=0;
    end
    if x(i)<0
        x(i)=0;
    end
    if Cw(i)<0
        Cw(i)=0;
    end
    Cw(i+1)=dt*(bw(i)*(1-Cw(i))*Cw(i)-B*x(i)*Cw(i))+Cw(i);
    x(i+1)=dt*(-F*x(i)+G*x(i)*(Cw(i+1)))+x(i);
    
end
% ta = (ta4)^(1/4)
subplot(2, 2, 1)
yyaxis left
plot(time, Cw, 'm')
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
legend('Fraction covered by white daisies', 'Herbivores per unit area')
%xlim([0 20])

subplot(2, 2, 2)
plot([time(1:end-1)], Tw, 'm')
title('Temperature')
xlabel('Time')
ylabel('Temp')
legend('White daisies Temp')

subplot(2, 2, 3)
plot(x, Cw, 'm')
xlabel('Herbivores per unit area')
ylabel('Fraction covered by white daisies')
hold on
% quiver(x,Cb,gradient(x),gradient(Cb))
title('Herbivore/white daisy phase space')