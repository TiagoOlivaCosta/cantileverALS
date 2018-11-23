function [ry,y] = inov(sys,x,z)
% Retorna a Inovacao a cada iteracao, relativa ao modelo e a medida
if ~isa(sys, 'ss') % a variavel sys deve ser um sistema no espaco de estados
  error('variable sys is not of state space class')
end

A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

% Inovacao
y = C*x;    % observacao
ry = z-y; % residuo da inovacao
end
