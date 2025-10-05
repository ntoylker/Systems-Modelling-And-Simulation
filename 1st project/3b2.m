clc; clear; close all;

% Πραγματικές παράμετροι συστήματος
m_true = 0.75;   % kg
L_true = 1.25;   % m
c_true = 0.15;   % N*m*sec
g = 9.81;        % m/sec^2
A0 = 4;          % Αρχικό πλάτος εισόδου
omega = 2;       % Συχνότητα εισόδου [rad/s]
T_sim = 20;      % Χρόνος προσομοίωσης (sec)

% Πίνακας περιόδων δειγματοληψίας
Ts_values = [0.0001 0.001 0.005 0.01 0.05 0.1 0.25 0.5];

% Αποθήκευση σφαλμάτων εκτίμησης
errors_m = zeros(size(Ts_values));
errors_L = zeros(size(Ts_values));
errors_c = zeros(size(Ts_values));

% Ορισμός συνεχούς συστήματος
A = [0 1; -g/L_true -c_true/(m_true*L_true^2)];
B = [0; 1/(m_true*L_true^2)];
C = [1 0];
D = 0;
sys = ss(A, B, C, D);

for i = 1:length(Ts_values)
    Ts = Ts_values(i); % Τρέχουσα περίοδος δειγματοληψίας

    % Δημιουργία χρόνου και εισόδου
    t_sampled = 0:Ts:T_sim;
    u_sampled = A0 * sin(omega * t_sampled);

    % Απόκριση συστήματος
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

    % Υπολογισμός σχετικών σφαλμάτων
    errors_m(i) = abs((m_hat - m_true) / m_true) * 100;
    errors_L(i) = abs((L_hat - L_true) / L_true) * 100;
    errors_c(i) = abs((c_hat - c_true) / c_true) * 100;
end

% Δημιουργία διαγραμμάτων για τα σφάλματα εκτίμησης των παραμέτρων
figure('Position', [100, 100, 900, 600]);

subplot(3,1,1);
semilogx(Ts_values, errors_m, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Περίοδος Δειγματοληψίας T_s [s]');
ylabel('Σφάλμα εκτίμησης m [%]');
title('Σφάλμα εκτίμησης της παραμέτρου m');

subplot(3,1,2);
semilogx(Ts_values, errors_L, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Περίοδος Δειγματοληψίας T_s [s]');
ylabel('Σφάλμα εκτίμησης L [%]');
title('Σφάλμα εκτίμησης της παραμέτρου L');

subplot(3,1,3);
semilogx(Ts_values, errors_c, 'go-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Περίοδος Δειγματοληψίας T_s [s]');
ylabel('Σφάλμα εκτίμησης c [%]');
title('Σφάλμα εκτίμησης της παραμέτρου c');
