clc; clear; close all;

% Ορισμός των παραμέτρων του συστήματος
m_true = 0.75;   % kg
L_true = 1.25;   % m
c_true = 0.15;   % N*m*sec
g = 9.81;        % m/sec^2
A0 = 4;          % Αρχικό πλάτος εισόδου
omega = 2;       % Συχνότητα εισόδου [rad/s]
T_sim = 20;      % Αυξημένος χρόνος προσομοίωσης (sec)
Ts = 0.1;        % Περίοδος δειγματοληψίας

% Διατύπωση του συνεχούς γραμμικού μοντέλου
A = [0 1; -g/L_true -c_true/(m_true*L_true^2)];
B = [0; 1/(m_true*L_true^2)];
C = [1 0];
D = 0;

sys = ss(A, B, C, D);

% Δημιουργία εισόδου και χρόνος
t_sampled = 0:Ts:T_sim;
u_sampled = A0 * sin(omega * t_sampled);

% Απόκριση με lsim
[y, ~, x] = lsim(sys, u_sampled, t_sampled);

% Φιλτράρισμα του y με Savitzky-Golay
y_smooth = sgolayfilt(y, 3, 11);

% Υπολογισμός 1ης και 2ης παραγώγου του y
y_dot = gradient(y_smooth, Ts);
y_ddot = gradient(y_dot, Ts);

% Μείωση δείγματος για αποφυγή θορύβου στα άκρα
N = length(y) - 5;

% Εφαρμογή μεθόδου ελαχίστων τετραγώνων στη μορφή της εξίσωσης κίνησης
Phi = [-y_smooth(1:N), -y_dot(1:N), u_sampled(1:N)'];
y_target = y_ddot(1:N);

% Εκτίμηση παραμέτρων
theta_hat = (Phi' * Phi + 1e-6 * eye(size(Phi,2))) \ (Phi' * y_target);
alpha1_hat = theta_hat(1);
alpha2_hat = theta_hat(2);
b1_hat = theta_hat(3);

% Υπολογισμός των L, m, c από τις εκτιμήσεις
L_hat = g / alpha1_hat;
m_hat = 1 / (b1_hat * L_hat^2);
c_hat = alpha2_hat * m_hat * L_hat^2;

% Προβολή αποτελεσμάτων εκτίμησης
fprintf('Εκτιμώμενες παράμετροι: \n');
fprintf('alpha1 = %.4f, alpha2 = %.4f, b1 = %.4f\n', alpha1_hat, alpha2_hat, b1_hat);
fprintf('L_hat = %.4f m, m_hat = %.4f kg, c_hat = %.4f N*m*sec\n', L_hat, m_hat, c_hat);

% Δημιουργία του εκτιμώμενου συστήματος
A_hat = [0 1; -alpha1_hat -alpha2_hat];
B_hat = [0; b1_hat];
C_hat = [1 0];
D_hat = 0;

sys_hat = ss(A_hat, B_hat, C_hat, D_hat);

% Υπολογισμός της εκτιμώμενης απόκρισης με lsim
y_est = lsim(sys_hat, u_sampled, t_sampled);

% Σύγκριση πραγματικής και εκτιμώμενης απόκρισης
figure;

subplot(2,1,1);
plot(t_sampled, y, 'b', 'LineWidth', 1.5);
hold on;
plot(t_sampled, y_est, 'r--', 'LineWidth', 1.5);
xlabel('Χρόνος (s)'); ylabel('\theta (rad)');
title('Σύγκριση πραγματικής και εκτιμώμενης γωνίας');
legend('Πραγματική', 'Εκτιμώμενη');
grid on;
hold on;
eq = y_smooth - y_est;

subplot(2,1,2);
plot(t_sampled, eq, 'k', 'LineWidth', 1.2);
xlabel('Χρόνος (s)');
ylabel('Σφάλμα eq(t) = q(t) - q̂(t)');
title('Σφάλμα εκτίμησης');
grid on;