function R = autoCovMat(y,i,D)
    % y -> vetor de medicao
    % i -> indice de medicao
    % D -> numero de delays

    N = length(y);
    if i>D
      for i=1:N
          Ci = y(:,end)*y(:,end-i)';
          %comment
      end
    end
