function [C, Ceq] = f_boot_nonlcon(x,LM,B,nbasis)
beta = x(1:nbasis);
phi = x(nbasis+1:end);
%phi = [0.5 0.5 1.1 0.2 0.1 0.6]';

% Make the smoothness contraint for compoennts
C = beta'*LM*beta+phi'*phi - B;
% And the equality constraints
Ceq = 0;




