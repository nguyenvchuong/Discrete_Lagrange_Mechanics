function ddth = singlePendulumDynamics(th,g,l)
%singlePendulumDynamics
%    DDTH = singlePendulumDynamics(TH,G,L)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    09-Jun-2023 11:50:58

ddth = -(g.*sin(th))./l;
