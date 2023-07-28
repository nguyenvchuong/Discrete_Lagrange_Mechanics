%% Get the result from EoM_Cart_Pole

syms x dx ddx th dth ddth 'real'; 
% x == horizontal position of the cart
% th == angle of the point-mass pendulum

syms m1 m2 k g l 'real';
% m1 == mass of the cart
% m2 == mass of the pendulum bob
% k = spring constant
% g == acceleration due to gravity
% l = length of the pendulum

%Generalized coordinates:
q = [x, th];
dq = [dx, dth];

M = [m1 + m2, l*m2*cos(th); 
    l*m2*cos(th),l^2*m2];

V =-g*l*m2*cos(th);

L = 0.5*dq*M*dq';

%% Derive discrete Lagrangian Mechanics
% syms qk0 qk qk1 % qk0=q(k-1), qk=q(k), qk1= q(k+1)
syms x0 th0 x1 th1 'real'
syms h % time step
qk0 = [x0, th0];
qk = [x, th];
qk1= [x1, th1];

M = [m1 + m2, l*m2*cos(th); 
    l*m2*cos(th),l^2*m2];

V =-g*l*m2*cos(th);

% L(q,q_dot) = 0.5*q_dot'*M*q_dot -V(q)

% Ld(qk,qk1) = 0.5*h*(L(qk,(qk1-qk)/h)+L(qk1, (qk1-qk)/h)
% L1_qk_qk1=L(qk,(qk1-qk)/h), L2_qk_qk1=L(qk1, (qk1-qk)/h)
L1_qk_qk1= 0.5/(h^2)*(qk1-qk)*M*(qk1-qk)'-V;
L2_qk_qk1= 0.5/(h^2)*(qk1-qk)*subs(M,th,th1)*(qk1-qk)'- subs(V,th,th1);
Ld_qk_qk1 = L1_qk_qk1 + L2_qk_qk1;
D1Ld_qk_qk1= jacobian(Ld_qk_qk1,qk);

% Ld(qk0,qk) = 0.5*h*(L(qk0,(qk-qk0)/h)+L(qk, (qk-qk0)/h)

L1_qk0_qk= 0.5/(h^2)*(qk-qk0)*subs(M,th,th0)*(qk-qk0)'-subs(V,th,th0);
L2_qk0_qk= 0.5/(h^2)*(qk-qk0)*M*(qk-qk0)'- V;
Ld_qk0_qk = L1_qk0_qk + L2_qk0_qk;
D2Ld_qk0_qk= jacobian(Ld_qk0_qk, qk);

ans = simplify(D1Ld_qk_qk1 + D2Ld_qk0_qk);










