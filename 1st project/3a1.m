% Καθαρισμός workspace
clc; clear; close all;

% Πραγματικές Παράμετροι Συστήματος
m_true = 0.75;
L_true = 1.25;
c_true = 0.15;
g = 9.81;
A0 = 4;
omega = 2;
Ts = 0.1;

% Προσομοίωση Συστήματος
A = [0 1; -g/L_true -c_true/(m_true*L_true^2)];
B = [0; 1/(m_true*L_true^2)];
odefun = @(t, x) A*x + B*A0*sin(omega*t);
tspan = 0:0.001:20;
x0 = [0; 0];
[t, x] = ode45(odefun, tspan, x0);

% Δειγματοληψία
t_sampled = (0:Ts:20)';
q_sampled = interp1(t, x(:,1), t_sampled, 'linear');
dq_sampled = interp1(t, x(:,2), t_sampled, 'linear');
u_sampled = A0 * sin(omega * t_sampled);

% Ορισμός Φίλτρων
lambda = 1;
s = tf('s');
filter_s_over_s1 = s / (s + lambda);
filter_1_over_s1 = 1 / (s + lambda);

% Διαφορετικά SNR
SNR_values = [10, 20, 30];
errors_m = zeros(size(SNR_values));
errors_L = zeros(size(SNR_values));
errors_c = zeros(size(SNR_values));

for i = 1:length(SNR_values)
    snr = SNR_values(i);
    
    % Προσθήκη θορύβου με awgn
    q_noisy = awgn(q_sampled, snr, 'measured');
    dq_noisy = awgn(dq_sampled, snr, 'measured');
    u_noisy = awgn(u_sampled, snr, 'measured');

    % Φιλτραρισμένες ποσότητες
    sy = gradient(q_noisy, Ts);
    s_over_s1_y = lsim(filter_s_over_s1, q_noisy, t_sampled);
    s2_over_s1_y = sy - s_over_s1_y;
    phi1 = -s_over_s1_y;
    phi2 = -lsim(filter_1_over_s1, q_noisy, t_sampled);
    phi3 = lsim(filter_1_over_s1, u_noisy, t_sampled);

    % Least Squares
    Y = s2_over_s1_y(10:end);
    Phi = [phi1(10:end), phi2(10:end), phi3(10:end)];
    theta = (Phi' * Phi + 1e-6 * eye(3)) \ (Phi' * Y);

    % Εκτίμηση παραμέτρων
    alpha1 = theta(1);
    alpha2 = theta(2);
    b1 = theta(3);

    L_est = g / alpha2;
    m_est = 1 / (b1 * L_est^2);
    c_est = alpha1 * m_est * L_est^2;

    % Υπολογισμός σφαλμάτων (%)
    errors_m(i) = abs((m_est - m_true) / m_true) * 100;
    errors_L(i) = abs((L_est - L_true) / L_true) * 100;
    errors_c(i) = abs((c_est - c_true) / c_true) * 100;
    
    % Εμφάνιση αποτελεσμάτων
    fprintf('\n--- SNR = %d dB ---\n', snr);
    fprintf('m εκτιμ.: %.4f (σφάλμα %.2f%%)\n', m_est, errors_m(i));
    fprintf('L εκτιμ.: %.4f (σφάλμα %.2f%%)\n', L_est, errors_L(i));
    fprintf('c εκτιμ.: %.4f (σφάλμα %.2f%%)\n', c_est, errors_c(i));
end

% Γράφημα Σφαλμάτων
figure('Position', [100 100 800 500]);
plot(SNR_values, errors_m, '-o', 'LineWidth', 2); hold on;
plot(SNR_values, errors_L, '-s', 'LineWidth', 2);
plot(SNR_values, errors_c, '-^', 'LineWidth', 2);
xlabel('SNR (dB)'); ylabel('Σφάλμα Εκτίμησης (%)');
title('Επίδραση Θορύβου στις Εκτιμήσεις Παραμέτρων');
legend('m', 'L', 'c');
grid on;


% Διάγραμμα Εκτιμώμενων Παραμέτρων vs Πραγματικές
figure('Name','Εκτιμώμενες Παράμετροι vs Πραγματικές','Position',[100 100 1000 400]);

% Δημιουργία πινάκων με τις εκτιμήσεις
m_est_all = zeros(size(SNR_values));
L_est_all = zeros(size(SNR_values));
c_est_all = zeros(size(SNR_values));

% Επαναυπολογισμός εκτιμήσεων για αποθήκευση (ή εναλλακτικά να τις αποθηκεύαμε στο αρχικό loop)
for i = 1:length(SNR_values)
    snr = SNR_values(i);
    
    q_noisy = awgn(q_sampled, snr, 'measured');
    dq_noisy = awgn(dq_sampled, snr, 'measured');
    u_noisy = awgn(u_sampled, snr, 'measured');

    sy = gradient(q_noisy, Ts);
    s_over_s1_y = lsim(filter_s_over_s1, q_noisy, t_sampled);
    s2_over_s1_y = sy - s_over_s1_y;
    phi1 = -s_over_s1_y;
    phi2 = -lsim(filter_1_over_s1, q_noisy, t_sampled);
    phi3 = lsim(filter_1_over_s1, u_noisy, t_sampled);

    Y = s2_over_s1_y(10:end);
    Phi = [phi1(10:end), phi2(10:end), phi3(10:end)];
    theta = (Phi' * Phi + 1e-6 * eye(3)) \ (Phi' * Y);

    alpha1 = theta(1);
    alpha2 = theta(2);
    b1 = theta(3);

    L_est = g / alpha2;
    m_est = 1 / (b1 * L_est^2);
    c_est = alpha1 * m_est * L_est^2;

    % Αποθήκευση
    m_est_all(i) = m_est;
    L_est_all(i) = L_est;
    c_est_all(i) = c_est;
end

% Υποπλοκές για κάθε παράμετρο
subplot(1,3,1);
plot(SNR_values, m_est_all, 'bo-', 'LineWidth', 2); hold on;
yline(m_true, 'r--', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)'); ylabel('m [kg]');
title('Εκτίμηση μάζας m');
legend('Εκτίμηση', 'Πραγματική');

subplot(1,3,2);
plot(SNR_values, L_est_all, 'gs-', 'LineWidth', 2); hold on;
yline(L_true, 'r--', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)'); ylabel('L [m]');
title('Εκτίμηση μήκους L');
legend('Εκτίμηση', 'Πραγματική');

subplot(1,3,3);
plot(SNR_values, c_est_all, 'm^-', 'LineWidth', 2); hold on;
yline(c_true, 'r--', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)'); ylabel('c [N·m·s]');
title('Εκτίμηση απόσβεσης c');
legend('Εκτίμηση', 'Πραγματική');

sgtitle('Σύγκριση Εκτιμώμενων και Πραγματικών Παραμέτρων');

% Σύγκριση Θορυβημένων vs Πραγματικών Δεδομένων για διάφορα SNR
SNR_values = [10, 20, 30]; % Δοκιμές για διαφορετικά επίπεδα θορύβου

figure('Name','Σύγκριση Θορυβημένων και Πραγματικών Δεδομένων','Position',[100, 100, 1200, 800]);

for i = 1:length(SNR_values)
    snr = SNR_values(i);
    
    % Προσθήκη θορύβου
    q_noisy = awgn(q_sampled, snr, 'measured');
    dq_noisy = awgn(dq_sampled, snr, 'measured');

    % Σύγκριση Πραγματικής vs Θορυβημένης απόκρισης q(t)
    subplot(3,2,2*i-1);
    plot(t_sampled(10:end), q_sampled(10:end), 'k', 'LineWidth', 1.5); hold on;
    plot(t_sampled(10:end), q_noisy(10:end), 'r--', 'LineWidth', 1.5);
    title(sprintf('Σύγκριση q(t) με SNR=%d dB', snr));
    legend('Πραγματική q(t)', 'Θορυβημένη q(t)');
    grid on;

    % Σύγκριση Πραγματικής vs Θορυβημένης ταχύτητας dq/dt
    subplot(3,2,2*i);
    plot(t_sampled(10:end), dq_sampled(10:end), 'k', 'LineWidth', 1.5); hold on;
    plot(t_sampled(10:end), dq_noisy(10:end), 'r--', 'LineWidth', 1.5);
    title(sprintf('Σύγκριση dq/dt με SNR=%d dB', snr));
    legend('Πραγματική dq/dt', 'Θορυβημένη dq/dt');
    grid on;
end
