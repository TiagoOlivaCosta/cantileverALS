% 0th order Takagi-Sugeno covariance estimation
% run setupFuzzy.m
 
 domainr = linspace(0,2);
 domainP = linspace(-2,2);
 mu1 = membr(domainr);
 mu2 = membP(domainP);

%Regras
           % r P u 
 ruleList = [1 3 3 1 1;...
             1 2 3 1 1;...
             1 1 3 1 1;...
            
             2 1 1 1 1;...
             2 3 5 1 1;...
             2 2 2 1 1;...
            ]; 
            
            
     
            
i1 = ruleList(:,1);
i2 = ruleList(:,2);
i3 = ruleList(:,3);

        
    out = [-8, -2, 0, 2, 8];
for i=1:length(mu1)
    for j = 1:length(mu2)
         antecedent = min(mu1(i1,i),mu2(i2,j));
        
              consequent = antecedent'*out(i3)';
            
        if  sum(antecedent)>0
          y = consequent/sum(antecedent);
        else
            y = 0;
        end

         
         z(j,i) = y;
     end
end
figure
mesh(domainr', domainP', z)
% mesh(z)
title('Acao de controle')
xlabel('|r|')
ylabel('(Pk-Pr)/Pr')
zlabel('u')

figure
mesh(domainr', domainP', abs(z))
% mesh(z)
title('Modulo da acao de controle')
xlabel('|r|')
ylabel('(Pk-Pr)/Pr')
zlabel('u')