% Extend Kalman Filter for Tricycle Mobile Robot Type 2
clc;
close all;

G = 13;
% Measurement data
measurementdata = load('measdat13.mat');
Z = measurementdata.me;

% Control signal 
controlsignal = load('contrlsig.mat');
V = controlsignal.inp;

d = sqrt(2.5 + 0.1*G);
dt = 0.05;% time step
t = 0: dt : 40; % total run time
control = [0.1; 0.2];

%Position 
p_pre = zeros(3, length(V) + 1);
p_actual = zeros(3, length(V) + 1);
p_kalman = zeros(3, length(V) + 1);
X_ = zeros(3, length(V) + 1);

% Compute 
%Actual
for i = 2: length(t)
    p_actual(1,i) = p_actual(1,i-1) + V(1,i-1)*cos(p_actual(3,i-1))*dt;
    p_actual(2,i) = p_actual(2,i-1) + V(1,i-1)*sin(p_actual(3,i-1))*dt;
    p_actual(3,i) = p_actual(3,i-1) + V(1,i-1)*tan(V(2,i-1))/d*dt;
end
%Prediction
for i = 2: length(t)
    p_pre(1,i) = p_pre(1, i-1) + control(1)*cos(p_pre(3, i-1))*dt;
    p_pre(2,i) = p_pre(2, i-1) + control(1)*sin(p_pre(3, i-1))*dt;
    p_pre(3,i) = p_pre(3, i-1) + control(1)*tan(control(2))/d*dt;
end

%Kalman filter
LM1 = [3; 4];
LM2 = [4;4];

Q = [(2.5*10^-3) 0; 0 (3.6*10^-3)];   % variance of the input noise.

R = [(10^-6) 0 0 0;...
    0 (10^-6) 0 0;...
    0 0 (7.62*10^-5) 0;...
    0 0 0 (7.62*10^-5)];  % variance of the measurement noise.

P = zeros(3);   

for i = 2: length(t)
    %Time update (Predict)
    A = [1 0 -control(1)*sin(p_kalman(3,i-1))*dt;...
        0 1 control(1)*cos(p_kalman(3,i-1))*dt;...
        0 0 1];
    
    W = [cos(p_kalman(3,i-1))*dt 0;...
         sin(p_kalman(3,i-1))*dt 0;...
         tan(control(2))*dt/d control(1)*dt/(d*(cos(control(2)))^2)];
     P_ = A*P*A' + W*Q*W';
     
     %Measurement update (Correct) 
     measurement = find_h_matrix([p_kalman(:,i-1); LM1 ; LM2]);
     
     H = [(p_kalman(1,i-1)-LM1(1))/measurement(1) (p_kalman(2,i-1)-LM1(2))/measurement(1) 0;...
       (p_kalman(1,i-1)-LM2(1))/measurement(2) (p_kalman(2,i-1)-LM2(2))/measurement(2) 0;...
       (LM1(2)-p_kalman(2,i-1))/measurement(1)^2 (LM1(1)-p_kalman(1,i-1))/measurement(1)^2 -1;...
       (LM2(2)-p_kalman(2,i-1))/measurement(1)^2 (LM2(1)-p_kalman(1,i-1))/measurement(1)^2 -1];

   %Kalman gain
     K = P_*H'*(H*P_*H' + R)^(-1);
     
     X_(1,i) = p_actual(1, i-1) + control(1)*cos(p_actual(3, i-1))*dt;
     X_(2,i) = p_actual(2, i-1) + control(1)*sin(p_actual(3, i-1))*dt;
     X_(3,i) = p_actual(3, i-1) + control(1)*tan(control(2))/d*dt;
     
     h = find_h_matrix([X_(:,i-1); LM1 ; LM2]);
     
     p_kalman(:,i) = X_(:, i) + K*(Z(:,i) - h);
     
    P = (eye(3) - K*H)*P_;
end 

% Plot result
figure(1);
hold on;
plot(p_kalman(1,:), p_kalman(2,:), '-k');
plot(p_pre(1,:), p_pre(2,:), '-r');
plot(p_actual(1,:), p_actual(2,:), '-b'); 
l = legend('Kalman','Dự đoán','Thực tế');
l.Location = 'southeast';
xlabel('X(m)')
ylabel('Y(m)');
title('�?ư�?ng đi thực và dự đoán của Robot');
grid on;

% Error
error_kalman = p_actual - p_kalman;
error_predict = p_actual - p_pre;

figure(2);
hold on;
plot(t, error_predict(1,:), '-k');
plot(t, error_kalman(1,:), '-b');
l = legend('Sai số dự đoán', 'Sai số Kalman');
l.Location = 'southwest';
xlabel('Th�?i gian(step)')
ylabel('Sai số(m)');
title('Sai số theo phương X');
grid on;

figure(3);
hold on;
plot(t, error_predict(2,:), '-k');
plot(t, error_kalman(2,:), '-b');
l = legend('Sai số dự đoán', 'Sai số Kalman');
l.Location = 'southwest';
xlabel('Th�?i gian(step)')
ylabel('sai số(m)');
title('Sai số theo phương Y');
grid on;

figure(4);
hold on;
plot(t, error_predict(1,:), '-k');
plot(t, error_kalman(1,:), '-b');
l = legend('Sai số dự đoán', 'Sai số Kalman');
l.Location = 'southwest';
xlabel('Th�?i gian(step)')
ylabel('Sai số(m)');
title('Sai số theo phương theta');
grid on;

