clear; clc; close all;

function dtheta = estimatorODE(t, theta, x_data, u)
    global dt;

    % Ανάκτηση παραμέτρων
    xhat = theta(1:2);
    Ahat = reshape(theta(3:6), 2, 2);
    Bhat = theta(7:8);

    % Παρεμβολή του πραγματικού x(t)
    xt = interp1(x_data.t, x_data.val', t, 'linear', 'extrap')';
    ut = u(t);

    % Σφάλμα
    e = xt - xhat;

    % Εκτιμητής κατάστασης (παράλληλη δομή)
    dxhat = Ahat * xhat + Bhat * ut;

    % Εκτίμηση παραμέτρων
    dAhat = xhat * e';     % 2x1 * 1x2 = 2x2
    dBhat = ut * e;        % scalar * 2x1 = 2x1

    % --- Ενσωμάτωση Προβολής ---
    % Προβλεπόμενες επόμενες τιμές
    A_next = Ahat + dt * dAhat;
    B_next = Bhat + dt * dBhat;

    % Προβολή περιορισμού a11 ∈ [-3, -1]
    if A_next(1,1) < -3
        dAhat(1,1) = (-3 - Ahat(1,1)) / dt;
    elseif A_next(1,1) > -1
        dAhat(1,1) = (-1 - Ahat(1,1)) / dt;
    end

    % Προβολή περιορισμού b2 ≥ 1
    if B_next(2) < 1
        dBhat(2) = (1 - Bhat(2)) / dt;
    end

    % Παράγωγος όλου του εκτιμητή
    dtheta = [dxhat;
              dAhat(:);
              dBhat];
end

% Πραγματικές τιμές του συστήματος
A_true = [-2.15  0.25;
          -0.75 -2.00];
B_true = [0; 1.5];

% Χρονικές παράμετροι
T_final = 40;
dt = 0.01;
time = 0:dt:T_final;
N = length(time);

% Συνάρτηση εισόδου
u = @(t) sin(5*t);

% Αρχική κατάσταση x
x0 = [0; 0];
% --- Προσομοίωση πραγματικού x(t) ---
ode_x = @(t, x) A_true * x + B_true * u(t);
[~, xsim] = ode45(ode_x, time, x0);
x = xsim';

% Δομή για μετάδοση προς ode εκτιμητών
x_data.t = time;
x_data.val = x;

% --- Αρχικές συνθήκες εκτιμήσεων ---
xhat0 = [1; -1];
Ahat0 = [-2, 0;
          0, 0];
Bhat0 = [0; 3];
theta0 = [xhat0; Ahat0(:); Bhat0];


% Global dt για χρήση στον εκτιμητή
global dt;
dt = 0.01;

% --- Εκτέλεση online estimation μέσω ode45 ---
[t_out, theta_out] = ode45(@(t, th) estimatorODE(t, th, x_data, u), time, theta0);

% --- Εξαγωγή αποτελεσμάτων ---
x_hat = zeros(2, N);
A_hat = zeros(2, 2, N);
B_hat = zeros(2, N);

for i = 1:N
    th = theta_out(i, :)';
    x_hat(:, i) = th(1:2);
    A_hat(:,:,i) = reshape(th(3:6), 2, 2);
    B_hat(:, i) = th(7:8);
end

% --- Γραφήματα ---% --- Διάγραμμα σφάλματος e = x - x_hat ---
e_all = x - x_hat;

figure;
subplot(2,1,1);
plot(time, e_all(1,:), 'b', 'LineWidth', 1.5); title('Σφάλμα e_1(t)');
xlabel('t'); ylabel('e_1');

subplot(2,1,2);
plot(time, e_all(2,:), 'r', 'LineWidth', 1.5); title('Σφάλμα e_2(t)');
xlabel('t'); ylabel('e_2');

% --- Εκτιμήσεις Α ---
figure;
plot(time, squeeze(A_hat(1,1,:)), 'LineWidth', 1.5); hold on;
yline(A_true(1,1), 'k', 'LineWidth', 1.5);
title('\alpha_{11}'); xlabel('t'); ylabel('\alpha_{11}');

figure;
plot(time, squeeze(A_hat(1,2,:)), 'LineWidth', 1.5); hold on;
yline(A_true(1,2), 'k', 'LineWidth', 1.5);
title('\alpha_{12}'); xlabel('t'); ylabel('\alpha_{12}');

figure;
plot(time, squeeze(A_hat(2,1,:)), 'LineWidth', 1.5); hold on;
yline(A_true(2,1), 'k', 'LineWidth', 1.5);
title('\alpha_{21}'); xlabel('t'); ylabel('\alpha_{21}');

figure;
plot(time, squeeze(A_hat(2,2,:)), 'LineWidth', 1.5); hold on;
yline(A_true(2,2), 'k', 'LineWidth', 1.5);
title('\alpha_{22}'); xlabel('t'); ylabel('\alpha_{22}');

% --- Εκτιμήσεις B ---
figure;
plot(time, B_hat(1,:), 'LineWidth', 1.5); hold on;
yline(B_true(1), 'k', 'LineWidth', 1.5);
title('\beta_1'); xlabel('t'); ylabel('\beta_1');

figure;
plot(time, B_hat(2,:), 'LineWidth', 1.5); hold on;
yline(B_true(2), 'k', 'LineWidth', 1.5);
title('\beta_2'); xlabel('t'); ylabel('\beta_2');

% --- Εκτιμήσεις Καταστάσεων ---
figure;
subplot(2,1,1);
plot(time, x(1,:), 'b', 'LineWidth', 1.5); hold on;
plot(time, x_hat(1,:), 'r', 'LineWidth', 1.5);
legend('x_1', '\hat{x}_1'); title('Κατάσταση x_1 και Εκτίμηση');

subplot(2,1,2);
plot(time, x(2,:), 'b', 'LineWidth', 1.5); hold on;
plot(time, x_hat(2,:), 'r', 'LineWidth', 1.5);
legend('x_2', '\hat{x}_2'); title('Κατάσταση x_2 και Εκτίμηση');
