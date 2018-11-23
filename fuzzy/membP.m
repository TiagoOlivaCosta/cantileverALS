function P = membP(domain)
close all
x = domain;
ZP = domain(end);
ZN = domain(1);


%ZE
si = (ZP)/3.75;
ci = (ZP+ZN)/2;
ze = gaussmf(x, [si, ci]);

ZP = gaussmf(x, [si, 1.1*si]);
ZN  = gaussmf(x, [si, -1.1*si]);

inc = 3.75/(si); 
ce = 1.15*si;
p = sigmf(x, [inc, 2*ce]);
n = sigmf(x, [-inc, -2*ce]);

comp = 1-ze;
% p = [zeros(1,length(domain)/2-1), comp(length(comp)/2:end)];
% n = [comp(1:length(comp)/2), zeros(1,length(domain)/2) ];

P = [n; ZN; ze; ZP; p];
% 
%  figure
%  plot(x, n)
% hold on
% plot(x, ze)
% plot(x, p)
% plot(x, ZN)
% plot(x, ZP)
% figure
% plot(x,sum(P))


end