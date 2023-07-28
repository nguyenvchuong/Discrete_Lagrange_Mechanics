%EoM_double link.m

% This script derives the equations of motion for a double pendulum (mimic 2 DoF leg of the robot), using the lagrange equations 
% and the matlab symbolic toolbox.

% The lagrangian (L) is defined as:
% 
%       L = T - U
%
% where 
%       T = system's kinetic energy
%       U = system's potential energy

% How we go about expressing the equations of motion:
%
%       DL      D  / DL \         * Note that some of those 'D' should be
%       ---  =  -- | -- |           curvy 'D' to represent partial
%       Dq      Dt \ Dq /           derivatives

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        set up variables                                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
clc; clear;

syms th1 dth1 ddth1 th2 dth2 ddth2 'real'; 
% th1 == thigh angle
% th2 == calf angle

syms m1 m2 Ic1 Ic2 g L1 L2 'real';
% m1 == mass of the thigh
% m2 == mass of the calf
% g == acceleration due to gravity
% l = length of the pendulum

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                            vector stuff                                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

p1 = [-0.5*L1*sin(th1); -0.5*L1*cos(th1)];
p2 = [-L1*sin(th1) - 0.5*L2*sin(th1 + th2); -L1*cos(th1)-0.5*L2*cos(th1 + th2)];

dp1 = [-0.5*L1*cos(th1)*dth1; 0.5*L1*sin(th1)*dth1];
dp2 = [-L1*cos(th1)*dth1 - 0.5*L2*cos(th1 + th2)*(dth1 + dth2); L1*sin(th1)*dth1 + 0.5*L2*sin(th1+th2)*(dth1 + dth2)];
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                          Lagrangian Definitions                         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%Kinetic energy:
T = (1/2)*m1*sum(dp1.*dp1) + (1/2)*m2*sum(dp2.*dp2) + 1/2 *(1/12*m1*L1^2)* dth1^2 + 1/2 *(1/12*m2*L2^2)* (dth1+dth2)^2;

%Potential energy:
U = m1*g*p1(2)+ m2*g*p2(2);

%Lagrangian:
L = T - U; 

%Generalized coordinates:
q = [th1, th2];
dq = [dth1, dth2];
ddq = [ddth1, ddth2];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                  evaluate partial derivatives                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%              DL  
%  DL_Dq  ==  ---      Note that 'D' is partial derivative here 
%              Dq
%
DL_Dq = jacobian(L,q')';

%              DL  
%  DL_Ddq  ==  ---      Note that 'D' is partial derivative here 
%              Ddq
%
DL_Ddq = jacobian(L,dq);

%                D  / DL  \         * Note that some of those 'd' should be
% DDL_DtDdq  ==  -- | --  |         curvy 'D' to represent partial
%                Dt \ Ddq /         derivatives
%
% Note the application of the chain rule:  (Quoting Andy Ruina: )
%      d BLAH / dt  =  jacobian(BLAH, [q qdot])*[qdot qddot]'
%
DDL_DtDdq = jacobian(DL_Ddq',[q, dq]) * [dq, ddq]';


%Write out as single equation and simplify:
EoM = -DL_Dq + DDL_DtDdq;
EoM = simplify(EoM);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                   mass matrix gymnastics                                %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% This is a confusing step - read carefully:
%
% We know that our equations are of the form:
%             EoM = M(q,dq)*qdd + f(q,dq) == 0;
%
% Thus, we can find M(q,dq) by using the jacobian command:

M = jacobian(EoM,ddq);


% Now, we want to find f(q,dq). We can do this by substituting in zero for
% the acceleration vector (dqq)

f = subs(EoM,ddq,sym([0, 0]));


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                               write files                               %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

writeDoubleLinkDynamics(f,M);

