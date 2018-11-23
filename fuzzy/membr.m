function r = membr(domain)
x = domain;
ZP = domain(end);
ZN = domain(1);

inc = 3*ZP; 
ce = (ZP-ZN)/2;
b = sigmf(x, [inc, ce]);
s = sigmf(x, [-inc, ce]);

%   plot(x, pn)
% hold on
% plot(x, pb)

r = [s; b];
    
end