function [energy,kinetic,potential] = singlePendulumEnergy(th,w,m,g,l)
%singlePendulumEnergy
%    [ENERGY,KINETIC,POTENTIAL] = singlePendulumEnergy(TH,W,M,G,L)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    09-Jun-2023 11:50:58

t2 = cos(th);
t3 = l.^2;
t4 = w.^2;
t5 = t2-1.0;
t7 = (m.*t3.*t4)./2.0;
t6 = g.*l.*m.*t5;
t8 = -t6;
energy = t7+t8;
if nargout > 1
    kinetic = t7;
end
if nargout > 2
    potential = t8;
end
