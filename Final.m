% S
syms s
m_1 = 2;
m_2 = 3.5;
r = 0.061;
I_1 = 0.32;
L = 0.4;
theta_2 = 0.01;
I_2 = 0.0065;
g = 9.81;
b = 0.001;
a = m_2*g*L;
D = 0;

H_1 = ((m_1+m_2)*(r^2)) + I_1;
H_2 = m_2*r*L*cos(theta_2);
H_3 = (m_2*(L^2))+I_2;

A = [-b*(H_1 + H_2 + (2*H_3)) / ((H_1*H_2) - (H_3^2)), (a*H_1) / ((H_1*H_2) - (H_3^2)), a / ((H_1*H_2) - (H_3^2));
         1, 0, 0;
         0, 1, 0];


B = [1;0;0];

C = [0, (H_1+H_3)/((H_1*H_2)-(H_3^2)),0];

P = [B,A*B,(A^2)*B];

disp(P);
rank(P);

Q = [C;
    C*A;
    C*(A^2)];

disp(Q);
rank(Q);
% Stability Analysis
% 1. Lyapunov Stability
eigenvalues_A = eig(A);

if all(real(eigenvalues_A) < 0)
    disp('The system is Lyapunov stable');
else
    disp('The system is not Lyapunov stable.');
end

% 2. BIBO Stability
D_mat = 0; 
transfer_function = ss(A, B, C, D);
[bibo_numerator, bibo_denominator] = tfdata(transfer_function);
if all(real(roots(bibo_denominator{1})) < 0)
    disp('The system is BIBO stable ');
else
    disp('The system is not BIBO stable');
end
sys = ss(A,B,C,D);
step(sys)
Tf = tf(sys);
rlocus(Tf)

desired_poles = [-2 + 2j, -2 - 2j, -5];
K = place(A, B, desired_poles);
A_cl = A - B * K;
sys_cl = ss(A_cl, B, C, D);
tf2 = tf(sys_cl);
rlocus(tf2)

figure;
step(sys_cl);
title('Closed-Loop Step Response');
grid on;

t = 0:0.01:10; 

x0 = [1; 0; 0]; 
A_cl = A - B * K; 
sys_cl = ss(A_cl, B, eye(3), 0);
[states, time] = initial(sys_cl, x0, t);

figure;
subplot(3, 1, 1);
plot(time, states(:, 1), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('State x_1');
title('State x_1 Response');
grid on;

subplot(3, 1, 2);
plot(time, states(:, 2), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('State x_2');
title('State x_2 Response');
grid on;

subplot(3, 1, 3);
plot(time, states(:, 3), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('State x_3');
title('State x_3 Response');
grid on;

sgtitle('Individual State Responses');