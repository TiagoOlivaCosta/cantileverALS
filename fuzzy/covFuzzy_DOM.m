function  [ R, Q] = covFuzzy(r, R, Q, Cp, S)
%Regras
           %  P u 
 ruleList = [ 1 5 ;...
              2 4 ;...            
              3 3 ;...
              4 2 ;...
              5 1 ;...
            ];          

i2 = ruleList(:,1);
i3 = ruleList(:,2);
out = [-5, -1, 0, 1, 5]/25; 

% calculo das funcoes de pertinencia
 S = diag(S);
 Cp = diag(Cp);
 DoM = S-Cp;

  L1 =  [[-5 5]*10^-11;
        [-5, 5]*10^-7];     
    domainP = linspace(L1(1,1), L1(1,2));
     memb2 = membP(domainP);
 for i = 1:length(DoM)

     m = round((DoM(i)-L1(i,1))*100/(L1(i,2)-L1(i,1)));
        
     if m>100
         m = 100;
     end
     if m < 1
         m = 1;
     end
   
     mu = memb2(:, m);
        
     % processamento fuzzy
         antecedent = min(1, mu(i2));  
         consequent = antecedent'*out(i3)';
            
        if  sum(antecedent)>0
          u = consequent/sum(antecedent);
        else
            u = 0;
        end

         
%         if abs(u)<0.01
%             u = 0;
%         end
        % atualizacao dos elementos da diagonal de R e Q
         Rn(i) = (1+u)*R(i,i);      
         %Qn(i) = (1-u)*Q(i,i);
        

 end 
  R = diag(Rn);
  %Q = diag([Qn, Qn]);
 
end
