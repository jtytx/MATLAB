function [C, Ceq] = f_nonlcon2(x,beta,LM,B)


% Make the smoothness constraint for the components in Theta
C = beta'*LM*beta+x'*x - B;

% And the equality constraints
Ceq = 0;