% Given parameters
m = 0.75;       % kg
L = 1.25;       % m
c = 0.15;       % N·m·s
g = 9.81;       % m/s^2
A0 = 4;         % N·m (amplitude of input torque)
omega = 2;      % rad/s (frequency of input torque)


% State-space matrices
A = [0 1;
     -g/L -c/(m*L^2)];

B = [0;
     1/(m*L^2)];

C = [1 0]; % Output: angle q(t)
D = 0;

% Time vector
dt = 1e-3; % Step size
tspan = 0:dt:20; % 20 seconds simulation

% Define input function u(t) = A0*sin(omega*t)
u = @(t) A0 * sin(omega * t);

% ODE function
ode_func = @(t, x) A*x + B*u(t);

% Initial conditions
x0 = [0; 0]; % q(0) = 0, dq/dt(0) = 0


% Solve ODE using an explicit function
[t, x] = ode45(@(t, x) pendulum_ode(t, x, A, B, A0, omega), tspan, x0);

% Extract outputs
q = x(:,1); % Extract angular displacement

% Plot results
figure;
subplot(2,1,1);
plot(t, q, 'b', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Angle q(t) [rad]');
title('Pendulum Angle Response');
grid on;

subplot(2,1,2);
plot(t, A0*sin(omega*t), 'r', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Control Input u(t) [N·m]');
title('Input Torque u(t)');
grid on;

% ODE function
function dxdt = pendulum_ode(t, x, A, B, A0, omega)
    u = A0 * sin(omega * t); % Define input torque
    dxdt = A*x + B*u; % Compute dx/dt
end