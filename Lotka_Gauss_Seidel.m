clear all
tot_time=100;
dt = 0.03;
Niter=ceil(tot_time/dt);

% Initial conditions for Lotka-Volterra model:
P = 1; % Growth rate for prey
B = 0.1; % Predation rate
F = 1.5; % Death rate for predator
G = 0.03; % Growth rate for predator
x(1) = 30;
y(1) = 3; 

% Forward Euler (unstable):
for i=1:Niter
      time(i+1) = dt*i;
     x(i+1) = dt*(P*x(i)-B*x(i)*y(i))+x(i);
     y(i+1) = dt*(-F*y(i)+G*x(i+1)*y(i))+y(i);
end

subplot(1, 2, 1)
plot(time, x, 'Linewidth', 2)
hold on
plot(time, y, 'Linewidth', 2)
% title('Populations of predator and prey per unit area')
xlabel('Time', 'Fontsize', 40)
ylabel('Number per unit area', 'Fontsize', 40)
legend('Number of Prey', 'Number of Predator', 'Fontsize', 20)
title('Prey and Predator over Time', 'Fontsize', 40)
xlim([0 50])
subplot(1, 2, 2)
plot(x, y, 'Linewidth', 2)
xlabel('Number of Prey', 'Fontsize', 40)
ylabel('Number of Predator', 'Fontsize', 40)
title('Prey/Predator Phase Space', 'Fontsize', 40)