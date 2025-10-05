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

% Προσομοίωση Συστήματος με ODE45
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
q_sampled  = interp1(t, x(:,1), t_sampled, 'linear', 'extrap');
dq_sampled = interp1(t, x(:,2), t_sampled, 'linear', 'extrap');
u_sampled = A0 * sin(omega * t_sampled);

% Προσθήκη Θορύβου στα Μετρήσιμα Σήματα (q, u)
SNR_values = [10, 20, 30];
for snr = SNR_values
    q_noisy = awgn(q_sampled, snr, 'measured');
    u_noisy = awgn(u_sampled, snr, 'measured');
    
    % Υπολογισμός Παραγώγου (Προσέγγιση dq/dt από q)
    dq_est = gradient(q_noisy, Ts);
    
    % Ορισμός Φίλτρων
    lambda = 1;  
    s = tf('s');
    filter_1_over_s1 = 1/(s + lambda);
    filter_s_over_s1 = s/(s + lambda);
    
    % Φιλτράρισμα Δεδομένων
    s_over_s1_y = lsim(filter_s_over_s1, q_noisy, t_sampled);
    phi1 = -s_over_s1_y;
    phi2 = -lsim(filter_1_over_s1, q_noisy, t_sampled);
    phi3 = lsim(filter_1_over_s1, u_noisy, t_sampled);
    
    % Διαμόρφωση προβλήματος ελαχίστων τετραγώνων
    Y = dq_est(10:end); % Χρησιμοποιούμε την εκτιμημένη ταχύτητα
    Phi = [phi1(10:end), phi2(10:end), phi3(10:end)];
    
    % Επίλυση για τις παραμέτρους
    theta = (Phi' * Phi + 1e-6 * eye(3)) \ (Phi' * Y);
    alpha1 = theta(1);  % = c/(m*L²)
    alpha2 = theta(2);  % = g/L
    b1 = theta(3);      % = 1/(m*L²)
    
    L_est = g / alpha2;
    m_est = 1 / (b1 * L_est^2);
    c_est = alpha1 * m_est * L_est^2;
    
    % Εμφάνιση αποτελεσμάτων
    fprintf('SNR = %d\n', snr);
    fprintf('m: Πραγματικό = %.3f kg, Εκτιμώμενο = %.3f kg\n', m_true, m_est);
    fprintf('L: Πραγματικό = %.3f m, Εκτιμώμενο = %.3f m\n', L_true, L_est);
    fprintf('c: Πραγματικό = %.3f N·m·s, Εκτιμώμενο = %.3f N·m·s\n', c_true, c_est);
    
    % Σύγκριση q_noisy με q_sampled
    figure;
    subplot(2,1,1);
    plot(t_sampled, q_sampled, 'b', 'LineWidth', 1.5); hold on;
    plot(t_sampled, q_noisy, 'r--', 'LineWidth', 1.5);
    title(['Σύγκριση q(t) με SNR = ', num2str(snr)]);
    legend('Πραγματικό q(t)', 'Θορυβημένο q(t)');
    grid on;
    
    subplot(2,1,2);
    plot(t_sampled, u_sampled, 'b', 'LineWidth', 1.5); hold on;
    plot(t_sampled, u_noisy, 'r--', 'LineWidth', 1.5);
    title(['Σύγκριση u(t) με SNR = ', num2str(snr)]);
    legend('Πραγματικό u(t)', 'Θορυβημένο u(t)');
    grid on;
end

