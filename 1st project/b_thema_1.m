% Καθαρισμός workspace
clc; clear; close all;

% Πραγματικές Παράμετροι Συστήματος
m_true = 0.75;      % Μάζα [kg]
L_true = 1.25;      % Μήκος [m]
c_true = 0.15;      % Συντελεστής απόσβεσης [N·m·s]
g = 9.81;           % Επιτάχυνση βαρύτητας [m/s²]
A0 = 4;             % Πλάτος εισόδου [N·m]
omega = 2;          % Συχνότητα εισόδου [rad/s]
Ts = 0.1;           % Περίοδος δειγματοληψίας [s]

% Προσομοίωση Συστήματος με ODE45 (Γραμμικοποιημένο μοντέλο)
A = [0 1; 
    -g/L_true -c_true/(m_true*L_true^2)];
B = [0; 
     1/(m_true*L_true^2)];
odefun = @(t, x) A*x + B*A0*sin(omega*t);
tspan = 0:0.001:20;  % Λεπτό βήμα για ακρίβεια
x0 = [0; 0];         % Αρχικές συνθήκες
[t, x] = ode45(odefun, tspan, x0);

% Δειγματοληψία σε Ts = 0.1 sec
t_sampled = (0:Ts:20)';
q_sampled = interp1(t, x(:,1), t_sampled, 'linear', 'extrap');
dq_sampled = interp1(t, x(:,2), t_sampled, 'linear', 'extrap');
u_sampled = A0 * sin(omega * t_sampled);

% Ορισμός Φίλτρων (με λ=1 σύμφωνα με το έγγραφο)
lambda = 1;  
s = tf('s');

% Φίλτρα για τους όρους:
% - s/(s+1)*y(t)
% - 1/(s+1)*y(t)
% - 1/(s+1)*u(t)
filter_s_over_s1 = s/(s + lambda);
filter_1_over_s1 = 1/(s + lambda);

% Υπολογισμός όρων με την τεχνική αποσύνθεσης
% 1. Υπολογισμός s*y (αριθμητική παράγωγος)
sy = gradient(q_sampled, Ts);

% 2. Υπολογισμός s/(s+1)*y (φιλτραρισμένη παράγωγος)
s_over_s1_y = lsim(filter_s_over_s1, q_sampled, t_sampled);

% 3. Τελικός όρος: s²/(s+1)*y = s*y - s/(s+1)*y
s2_over_s1_y = sy - s_over_s1_y;

% Υπολογισμός των φίλτρων για την παλινδρόμηση
phi1 = -s_over_s1_y;                     % -s/(s+1)*y
phi2 = -lsim(filter_1_over_s1, q_sampled, t_sampled);  % -1/(s+1)*y
phi3 = lsim(filter_1_over_s1, u_sampled, t_sampled);   % 1/(s+1)*u

% Διαμόρφωση προβλήματος ελαχίστων τετραγώνων
% Αφαίρεση μεταβατικού (πρώτα 10 δείγματα)
Y = s2_over_s1_y(10:end);
Phi = [phi1(10:end), phi2(10:end), phi3(10:end)];

% Επίλυση με κανονικοποίηση
theta = (Phi' * Phi + 1e-6 * eye(3)) \ (Phi' * Y);

% Ανάκτηση Φυσικών Παραμέτρων
alpha1 = theta(1);  % = c/(m*L²)
alpha2 = theta(2);  % = g/L
b1 = theta(3);      % = 1/(m*L²)

L_est = g / alpha2;
m_est = 1 / (b1 * L_est^2);
c_est = alpha1 * m_est * L_est^2;

% Εκτύπωση Αποτελεσμάτων
fprintf('=== Πραγματικές vs Εκτιμώμενες Παράμετροι ===\n');
fprintf('Πραγματικό m: %.3f kg, Εκτιμώμενο m: %.3f kg\n', m_true, m_est);
fprintf('Πραγματικό L: %.3f m, Εκτιμώμενο L: %.3f m\n', L_true, L_est);
fprintf('Πραγματικό c: %.3f N·m·s, Εκτιμώμενο c: %.3f N·m·s\n', c_true, c_est);

% Σχεδίαση Αποτελεσμάτων
figure('Position', [100, 100, 900, 600]);

% Σύγκριση q(t)
subplot(3,1,1);
plot(t_sampled(10:end), q_sampled(10:end), 'b', 'LineWidth', 1.5); hold on;
plot(t_sampled(10:end), lsim(ss(A,B,[1 0],0), u_sampled(10:end), t_sampled(10:end)), 'r--', 'LineWidth', 1.5);
title('Σύγκριση Πραγματικής και Εκτιμώμενης Απόκρισης q(t)');
legend('Πραγματική q(t)', 'Εκτιμώμενη q(t)');
grid on;

% Σύγκριση dq/dt
subplot(3,1,2);
plot(t_sampled(10:end), dq_sampled(10:end), 'b', 'LineWidth', 1.5); hold on;
plot(t_sampled(10:end), lsim(ss(A,B,[0 1],0), u_sampled(10:end), t_sampled(10:end)), 'r--', 'LineWidth', 1.5);
title('Σύγκριση Πραγματικής και Εκτιμώμενης Ταχύτητας dq/dt');
legend('Πραγματική dq/dt', 'Εκτιμώμενη dq/dt');
grid on;

% Σφάλμα εκτίμησης
subplot(3,1,3);
error = q_sampled(10:end) - lsim(ss(A,B,[1 0],0), u_sampled(10:end), t_sampled(10:end));
plot(t_sampled(10:end), error, 'k', 'LineWidth', 1.5);
title('Σφάλμα Εκτίμησης q(t)');
xlabel('Χρόνος [s]');
grid on;