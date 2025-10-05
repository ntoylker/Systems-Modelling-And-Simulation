clear; clc; close all;

% --- Πραγματικό Σύστημα ---
A_true = [-2.15 0.25; -0.75 -2.00];
B_true = [0; 1.5];

% --- Χρονικά ---
T_final = 20;
dt = 0.01;
time = 0:dt:T_final;
N = length(time);

% --- Είσοδος & διαταραχή ---
u = @(t) sin(5*t);
w_max = 1;
omega = @(t) rand() * w_max * rand_unit_vec();  % ||omega(t)|| ∈ [0, w_max]


% --- Προσομοίωση x(t) ---
x0 = [0; 0];
ode_x = @(t, x) A_true * x + B_true * u(t) + omega(t);
[~, xsim] = ode45(ode_x, time, x0);
x = xsim';
x_data.t = time;
x_data.val = x;

% --- Αρχικές εκτιμήσεις ---
xhat0 = [1; -1];
Ahat0 = [-2.0, 0; 0, 0];
Bhat0 = [0; 3.0];
theta0 = [xhat0; Ahat0(:); Bhat0];

% --- Παράμετροι μάθησης ---
g1 = 5; g2 = 5; gB = 5;
sigma = 1;

% --- Εκτιμητής (παράλληλη δομή) ---
estimatorODE = @(t, theta) estimator_parallel(t, theta, x_data, u, dt, ...
    g1, g2, gB, sigma);

[t_out, theta_out] = ode45(estimatorODE, time, theta0);

% --- Αποθήκευση αποτελεσμάτων ---
x_hat = zeros(2, N);
A_hat = zeros(2, 2, N);
B_hat = zeros(2, N);

for i = 1:N
    th = theta_out(i,:)';
    x_hat(:,i) = th(1:2);
    A_hat(:,:,i) = reshape(th(3:6), 2, 2);
    B_hat(:,i) = th(7:8);
end

% --- Υπολογισμός και αποθήκευση ω(t) σε κάθε t ---
omega_all = zeros(2, N);
for k = 1:N
    omega_all(:,k) = omega(time(k));
end

norm_omega = vecnorm(omega_all);  % ||ω(t)|| για κάθε t

% --- Διάγραμμα ||ω(t)|| και ω_max ---
figure;
plot(time, norm_omega, 'm', 'LineWidth', 2); hold on;
yline(w_max, 'k--', 'LineWidth', 2);
xlabel('t'); ylabel('||\omega(t)||');
title('Νόρμα Διαταραχής Πόλωσης \omega(t)');
legend('||\omega(t)||', '\omega_{max}');
grid on;


e = x - x_hat;

% --- Γραφήματα ---%
% --- Νορμοποιημένο σφάλμα ||e(t)|| ---
norm_e = vecnorm(e);  % Υπολογίζει τη νόρμα κάθε στήλης του e (δηλ. για κάθε t)

figure;
plot(time, norm_e, 'k', 'LineWidth', 2);
title('||e(t)|| - Νόρμα Σφάλματος Εκτίμησης');
xlabel('t'); ylabel('||e(t)||');
grid on;

figure;
subplot(2,1,1); plot(time, e(1,:), 'b', 'LineWidth', 2); title('e_1(t)');
subplot(2,1,2); plot(time, e(2,:), 'r', 'LineWidth', 2); title('e_2(t)');

figure; plot(time, squeeze(A_hat(1,1,:)), 'LineWidth', 2); hold on;
yline(A_true(1,1), 'k', 'LineWidth', 1.5); title('\alpha_{11}');

figure; plot(time, squeeze(A_hat(1,2,:)), 'LineWidth', 2); hold on;
yline(A_true(1,2), 'k', 'LineWidth', 1.5); title('\alpha_{12}');

figure; plot(time, squeeze(A_hat(2,1,:)), 'LineWidth', 2); hold on;
yline(A_true(2,1), 'k', 'LineWidth', 1.5); title('\alpha_{21}');

figure; plot(time, squeeze(A_hat(2,2,:)), 'LineWidth', 2); hold on;
yline(A_true(2,2), 'k', 'LineWidth', 1.5); title('\alpha_{22}');

figure; plot(time, B_hat(1,:), 'LineWidth', 2); hold on;
yline(B_true(1), 'k', 'LineWidth', 1.5); title('\beta_1');

figure; plot(time, B_hat(2,:), 'LineWidth', 2); hold on;
yline(B_true(2), 'k', 'LineWidth', 1.5); title('\beta_2');

figure;
subplot(2,1,1); plot(time, x(1,:), 'b', 'LineWidth', 2); hold on;
plot(time, x_hat(1,:), 'r', 'LineWidth', 2); title('x_1 και \hat{x}_1');
legend('x_1','\hat{x}_1');

subplot(2,1,2); plot(time, x(2,:), 'b', 'LineWidth', 2); hold on;
plot(time, x_hat(2,:), 'r', 'LineWidth', 2); title('x_2 και \hat{x}_2');
legend('x_2','\hat{x}_2');

% --- Συνάρτηση Εκτιμητή: Παράλληλη Δομή ---
function dtheta = estimator_parallel(t, theta, x_data, u, dt, g1, g2, gB, sigma)
    xhat = theta(1:2);
    Ahat = reshape(theta(3:6), 2, 2);
    Bhat = theta(7:8);

    xt = interp1(x_data.t, x_data.val', t, 'linear', 'extrap')';
    ut = u(t);
    e = xt - xhat;

    % Παράλληλη δυναμική
    dxhat = Ahat * xhat + Bhat * ut;

    % Εκτίμηση Α (ανά στοιχείο)
    da11 = g1 * e(1) * ut - g1 * sigma * Ahat(1,1);
    da12 = g1 * e(1) * ut - g1 * sigma * Ahat(1,2);
    da21 = g2 * e(2) * ut - g2 * sigma * Ahat(2,1);
    da22 = g2 * e(2) * ut - g2 * sigma * Ahat(2,2);
    dAhat = [da11 da12; da21 da22];

    % Εκτίμηση B με σ-τροποποίηση
    dBhat = gB * ut * e - sigma * Bhat;

    % --- Προβολές ---
    a11_next = Ahat(1,1) + dt * dAhat(1,1);
    if a11_next < -3.0
        dAhat(1,1) = (-3.0 - Ahat(1,1)) / dt;
    elseif a11_next > -1.0
        dAhat(1,1) = (-1.0 - Ahat(1,1)) / dt;
    end

    b2_next = Bhat(2) + dt * dBhat(2);
    if b2_next < 1.0
        dBhat(2) = (1 - Bhat(2)) / dt;
    end

    dtheta = [dxhat;
              dAhat(:);
              dBhat];
end

function w = rand_unit_vec()
    v = randn(2,1);
    w = v / norm(v + 1e-8);  % μοναδιαίο διάνυσμα
end
