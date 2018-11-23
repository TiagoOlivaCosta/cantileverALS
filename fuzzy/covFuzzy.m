function  [ R, Q, M] = covFuzzy(r, R, Q, Cp, S)
%Regras
           % r P u 
 ruleList = [1 3 3 1 1;...
             1 2 3 1 1;...
             1 1 3 1 1;...
            
             2 1 1 1 1;...
             2 3 5 1 1;...
             2 2 3 1 1;...
            ];          
i1 = ruleList(:,1);
i2 = ruleList(:,2);
i3 = ruleList(:,3);
out = [-5, -1, 0, 1, 5]/5; 

% calculo das funcoes de pertinencia
 DoM = S-Cp;
 r = abs(r);
  dominio = [0, 0.01, -0.1, 0.1]*10^-3;   
 r1 = dominio(1);
 r2 = dominio(2);
 P1 = dominio(3);
 P2 = dominio(4);
 
 domainr = linspace(0, 2);
 domainP = linspace(-1, 1);
 memb1 = membr(domainr);
 memb2 = membP(domainP);
 
%   if ~(size(R)==size(Q))
%              error('Erro. P e Q devem ter a mesma dimensao')
%   end
         
 for i = 1:length(DoM)
     m1 = round((r(i)-r1)*100/(r2-r1));
%      m1 = 100;
     m2 = round((DoM(i)-P1)*100/(P2-P1));
     
     if m1 >100
        m1 = 100;
     end
     if m1 < 1
             m1 = 1;
     end
     if m2>100
         m2 = 100;
     end
     if m2 < 1
         m2 = 1;
     end
     M = [m1, m2]';
     mu1 = memb1(:, m1);
     mu2 = memb2(:, m2);
        
     % processamento fuzzy
         antecedent = min(mu1(i1), mu2(i2));  
         consequent = antecedent'*out(i3)';
            
        if  sum(antecedent)>0
          u = consequent/sum(antecedent);
        else
            u = 0;
        end

         
   
        % atualizacao dos elementos da diagonal de R e Q
         Rn(i) = (1+u)*R(i,i);      
         Qn(i) = (1-u)*Q(i,i);
        

 end 
  R = diag(Rn);
  %Q = diag(Qn);
 
end
