function [xplus] = inov(sys,x,u)
% Retorna a Inovacao a cada iteracao, relativa ao modelo e a medida
if ~isa(sys, 'ss') % a variavel sys deve ser um sistema no espaco de estados
  error('variable sys is not of state space class')
end

  A = sys.A;
  B = sys.B;
  C = sys.C;
  D = sys.D;

xplus = A*x + B*u; % state propagation


end
