function [C, Ceq] = f_nonlcon1(x,phi,LM,B)


% Make the smoothness constraint for the components in Theta
C = x'*LM*x+phi'*phi - B;

% And the equality constraints
Ceq = 0;