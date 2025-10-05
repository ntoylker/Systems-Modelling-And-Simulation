clc; clear;

% Ορισμός των πραγματικών παραμέτρων του συστήματος
m_true = 0.75;   % kg
L_true = 1.25;   % m
c_true = 0.15;   % N*m*sec
g = 9.81;        % m/sec^2
omega = 2;       % Συχνότητα εισόδου [rad/s]
T_sim = 20;      % Χρόνος προσομοίωσης (sec)
Ts = 0.1;        % Περίοδος δειγματοληψίας

% Διατύπωση του συνεχούς γραμμικού μοντέλου
A = [0 1; -g/L_true -c_true/(m_true*L_true^2)];
B = [0; 1/(m_true*L_true^2)];
C = [1 0];
D = 0;

sys = ss(A, B, C, D);

% Επιλογή τιμών του A0 σε λογαριθμική κλίμακα
A0_values = logspace(-2, 2, 20);  % 20 σημεία στο διάστημα [0.01, 100]

% Αποθήκευση σφαλμάτων εκτίμησης
error_L = zeros(size(A0_values));
error_m = zeros(size(A0_values));
error_c = zeros(size(A0_values));

% Βρόχος για διαφορετικά A0
for i = 1:length(A0_values)
    A0 = A0_values(i);
    
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

    % Υπολογισμός των σφαλμάτων εκτίμησης
    error_L(i) = abs(L_hat - L_true) / L_true * 100;  % Σχετικό σφάλμα (%)
    error_m(i) = abs(m_hat - m_true) / m_true * 100;  % Σχετικό σφάλμα (%)
    error_c(i) = abs(c_hat - c_true) / c_true * 100;  % Σχετικό σφάλμα (%)
end

% Σχεδίαση σφαλμάτων εκτίμησης συναρτήσει του A0
figure;
semilogx(A0_values, error_L, 'b-o', 'LineWidth', 1.5); hold on;
semilogx(A0_values, error_m, 'r-s', 'LineWidth', 1.5);
semilogx(A0_values, error_c, 'g-d', 'LineWidth', 1.5);
xlabel('Πλάτος εισόδου A0');
ylabel('Σφάλμα εκτίμησης (%)');
title('Σφάλμα εκτίμησης συναρτήσει του A0');
legend('Σφάλμα L', 'Σφάλμα m', 'Σφάλμα c', 'Location', 'Best');
grid on;
