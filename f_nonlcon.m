function [C, Ceq] = f_nonlcon(x,LM,B,nbasis)

beta = x(1:nbasis);
%phi = x(nbasis+1:end);
phi = [0.1 0.6]';

% Make the smoothness constraint for the components in Theta
C = beta'*LM*beta+phi'*phi - B;

% And the equality constraints
Ceq = 0;