clear all
tot_time=100;
dt = 0.03;
Niter=ceil(tot_time/dt);

% Initial conditions for Lotka-Volterra model:
P = 1; % Growth rate for prey
B = 0.1; % Predation rate
F = 1.5; % Death rate for predator
G = 0.03; % Growth rate for predator
x(1) = 30; % 2 prey individuals per unit area
y(1) = 3; % 1 predator per unit area
% Z = 0.3 % extra prey death rate

% Forward Euler (unstable):
for i=1:Niter
      time(i+1) = dt*i;
     x(i+1) = dt*(P*x(i)-B*x(i)*y(i))+x(i);
     y(i+1) = dt*(-F*y(i)+G*x(i+1)*y(i))+y(i);
end

subplot(1, 2, 1)
plot(time, x)
hold on
plot(time, y)
% title('Populations of predator and prey per unit area')
xlabel('Time')
ylabel('Number per unit area')
legend('Number of Prey', 'Number of Predator')
title('Prey and Predator over Time')
xlim([0 50])
subplot(1, 2, 2)
plot(x, y)
xlabel('Number of Prey')
ylabel('Number of Predator')
title('Prey/Predator Phase Space')