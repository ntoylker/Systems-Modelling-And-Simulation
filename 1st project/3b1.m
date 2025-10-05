% Καθαρισμός workspace
clc; clear; close all;

% Πραγματικές Παράμετροι Συστήματος
m_true = 0.75;      % Μάζα [kg]
L_true = 1.25;      % Μήκος [m]
c_true = 0.15;      % Συντελεστής απόσβεσης [N·m·s]
g = 9.81;           % Επιτάχυνση βαρύτητας [m/s²]
A0 = 4;             % Πλάτος εισόδου [N·m]
omega = 2;          % Συχνότητα εισόδου [rad/s]
T_max = 20;         % Χρόνος προσομοίωσης [s]

% Πίνακας με διαφορετικές περιόδους δειγματοληψίας
Ts_values = [0.0001 0.001 0.005 0.01 0.05 0.1 0.25 0.5];

% Αποθήκευση σφαλμάτων εκτίμησης
errors_m = zeros(size(Ts_values));
errors_L = zeros(size(Ts_values));
errors_c = zeros(size(Ts_values));

% Προσομοίωση Συστήματος με ODE45
A = [0 1; 
    -g/L_true -c_true/(m_true*L_true^2)];
B = [0; 
     1/(m_true*L_true^2)];
odefun = @(t, x) A*x + B*A0*sin(omega*t);
tspan = 0:0.001:T_max;  % Λεπτό βήμα για ακρίβεια
x0 = [0; 0];            % Αρχικές συνθήκες
[t, x] = ode45(odefun, tspan, x0);

% Ορισμός Φίλτρων (λ=1)
lambda = 1;  
s = tf('s');
filter_s_over_s1 = s/(s + lambda);
filter_1_over_s1 = 1/(s + lambda);

% Βρόχος για διαφορετικές Ts
for i = 1:length(Ts_values)
    Ts = Ts_values(i);

    % Δειγματοληψία
    t_sampled = (0:Ts:T_max)';
    q_sampled = interp1(t, x(:,1), t_sampled, 'linear', 'extrap');
    dq_sampled = interp1(t, x(:,2), t_sampled, 'linear', 'extrap');
    u_sampled = A0 * sin(omega * t_sampled);

    % Υπολογισμός s²/(s+1)*y
    sy = gradient(q_sampled, Ts);
    s_over_s1_y = lsim(filter_s_over_s1, q_sampled, t_sampled);
    s2_over_s1_y = sy - s_over_s1_y;

    % Υπολογισμός όρων για παλινδρόμηση
    phi1 = -s_over_s1_y;
    phi2 = -lsim(filter_1_over_s1, q_sampled, t_sampled);
    phi3 = lsim(filter_1_over_s1, u_sampled, t_sampled);

    % Διαμόρφωση προβλήματος ελαχίστων τετραγώνων
    Y = s2_over_s1_y(10:end);
    Phi = [phi1(10:end), phi2(10:end), phi3(10:end)];

    % Επίλυση με κανονικοποίηση
    theta = (Phi' * Phi + 1e-6 * eye(3)) \ (Phi' * Y);

    % Ανάκτηση Φυσικών Παραμέτρων
    alpha1 = theta(1);  
    alpha2 = theta(2);  
    b1 = theta(3);      

    L_est = g / alpha2;
    m_est = 1 / (b1 * L_est^2);
    c_est = alpha1 * m_est * L_est^2;

    % Υπολογισμός σχετικών σφαλμάτων
    errors_m(i) = abs((m_est - m_true) / m_true) * 100;
    errors_L(i) = abs((L_est - L_true) / L_true) * 100;
    errors_c(i) = abs((c_est - c_true) / c_true) * 100;
end

% Σχεδίαση Σφαλμάτων Εκτίμησης Παραμέτρων
figure('Position', [100, 100, 900, 600]);

subplot(3,1,1);
semilogx(Ts_values, errors_m, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Περίοδος Δειγματοληψίας T_s [s]');
ylabel('Σφάλμα εκτίμησης μάζας [%]');
title('Σφάλμα εκτίμησης της παραμέτρου m');

subplot(3,1,2);
semilogx(Ts_values, errors_L, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Περίοδος Δειγματοληψίας T_s [s]');
ylabel('Σφάλμα εκτίμησης μήκους [%]');
title('Σφάλμα εκτίμησης της παραμέτρου L');

subplot(3,1,3);
semilogx(Ts_values, errors_c, 'go-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Περίοδος Δειγματοληψίας T_s [s]');
ylabel('Σφάλμα εκτίμησης απόσβεσης [%]');
title('Σφάλμα εκτίμησης της παραμέτρου c');
