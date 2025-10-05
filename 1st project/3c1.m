clc; clear; close all;

% Πραγματικές Παράμετροι Συστήματος
m_true = 0.75;      % kg
L_true = 1.25;      % m
c_true = 0.15;      % N·m·s
g = 9.81;           % m/s²
omega = 2;          % rad/s
Ts = 0.1;           % Περίοδος δειγματοληψίας [s]

% Πίνακας τιμών πλάτους A0
A0_values = logspace(-2, 2, 20); % 20 σημεία από 0.01 έως 100

% Αποθήκευση σφαλμάτων εκτίμησης
errors_m = zeros(size(A0_values));
errors_L = zeros(size(A0_values));
errors_c = zeros(size(A0_values));

% Ορισμός φίλτρων
lambda = 1;  
s = tf('s');
filter_s_over_s1 = s/(s + lambda);
filter_1_over_s1 = 1/(s + lambda);

% Προσομοίωση για κάθε τιμή του A0
for i = 1:length(A0_values)
    A0 = A0_values(i);
    
    % Ορισμός Συνεχούς Συστήματος
    A = [0 1; -g/L_true -c_true/(m_true*L_true^2)];
    B = [0; 1/(m_true*L_true^2)];
    odefun = @(t, x) A*x + B*A0*sin(omega*t);
    tspan = 0:0.001:20;  % Λεπτό βήμα για ακρίβεια
    x0 = [0; 0];         % Αρχικές συνθήκες
    [t, x] = ode45(odefun, tspan, x0);

    % Δειγματοληψία
    t_sampled = (0:Ts:20)';
    q_sampled = interp1(t, x(:,1), t_sampled, 'linear', 'extrap');
    dq_sampled = interp1(t, x(:,2), t_sampled, 'linear', 'extrap');
    u_sampled = A0 * sin(omega * t_sampled);

    % Υπολογισμός s*y και φιλτραρισμένων σημάτων
    sy = gradient(q_sampled, Ts);
    s_over_s1_y = lsim(filter_s_over_s1, q_sampled, t_sampled);
    s2_over_s1_y = sy - s_over_s1_y;

    % Υπολογισμός των φίλτρων για την παλινδρόμηση
    phi1 = -s_over_s1_y;
    phi2 = -lsim(filter_1_over_s1, q_sampled, t_sampled);
    phi3 = lsim(filter_1_over_s1, u_sampled, t_sampled);

    % Αφαίρεση μεταβατικών φαινομένων
    Y = s2_over_s1_y(10:end);
    Phi = [phi1(10:end), phi2(10:end), phi3(10:end)];

    % Επίλυση με ελαχιστοποίηση τετραγώνων
    theta = (Phi' * Phi + 1e-6 * eye(3)) \ (Phi' * Y);

    % Ανάκτηση φυσικών παραμέτρων
    alpha1 = theta(1);  
    alpha2 = theta(2);  
    b1 = theta(3);

    L_est = g / alpha2;
    m_est = 1 / (b1 * L_est^2);
    c_est = alpha1 * m_est * L_est^2;

    % Υπολογισμός σχετικών σφαλμάτων (%)
    errors_m(i) = abs((m_est - m_true) / m_true) * 100;
    errors_L(i) = abs((L_est - L_true) / L_true) * 100;
    errors_c(i) = abs((c_est - c_true) / c_true) * 100;
end

figure;
semilogx(A0_values, errors_L, 'b-o', 'LineWidth', 1.5); hold on;
semilogx(A0_values, errors_m, 'r-s', 'LineWidth', 1.5);
semilogx(A0_values, errors_c, 'g-d', 'LineWidth', 1.5);
xlabel('Πλάτος εισόδου A0');
ylabel('Σφάλμα εκτίμησης (%)');
title('Σφάλμα εκτίμησης συναρτήσει του A0');
legend('Σφάλμα L', 'Σφάλμα m', 'Σφάλμα c', 'Location', 'Best');
grid on;
