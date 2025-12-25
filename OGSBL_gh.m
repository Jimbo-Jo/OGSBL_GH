function [Pm, search_area] = OGSBL_gh(Y,dd,lambda,search_area,etc)

%initialization
[M,T]=size(Y);
A_theta = exp(-1i*2*pi*(0:M-1)'*(dd/lambda)*cosd(search_area));
N = length(search_area);
beta=(M*T)/(0.1*norm(Y, 'fro')^2);
gamma=mean(abs(A_theta'*Y),2);
converged = false;
iter = 0;

%parameters
maxiter = 500;
tol = 1e-3;
b = 1e-6;
scale_lambda = -2;
k_a1 = -scale_lambda/2 + 1;
k_a2 = 1e-6;
gamma_b0 = b*ones(N,1);
gamma_lambda0 = scale_lambda*ones(N,1);
c = 1e-6;
d = 1e-6;

while ~converged
      
      iter = iter + 1;
      gamma_last = gamma;
      
      %update mu and Sigma
      sigmaY = eye(M)/beta + A_theta*diag(gamma)*A_theta';
      mu = diag(gamma)*A_theta'/sigmaY*Y;
      sigma = diag(gamma) - diag(gamma)*A_theta'/sigmaY*A_theta*diag(gamma);
      
      %update gamma
      gamma_b = gamma_b0 + 2*diag(mu*mu'+ T*sigma);
      gamma_lambda = gamma_lambda0 - T;
      gamma_a = (k_a1+gamma_lambda0/2)./(k_a2 + gamma/2);
      gamma = sqrt(gamma_b./gamma_a).* besselk(gamma_lambda+1,sqrt(gamma_a.*gamma_b))./besselk(gamma_lambda,sqrt(gamma_a.*gamma_b));
      
      
      %update beta
      resid = Y-A_theta*mu;
      beta = (M*T+c)/(norm(resid, 'fro')^2+T*trace(A_theta*sigma*A_theta')+d);
      
      %stopping criteria
      erro=norm(gamma - gamma_last)/norm(gamma_last);
      if erro < tol || iter >= maxiter
          converged = true;
      end
      
end
      
      % offgrid DOA search
      [~, index_temp] = sort(gamma, 'descend');
      index = index_temp(1:etc);
      
      for j = 1:length(index)
          ii = index(j);   
          if ii ~= 1&&ii ~=length(search_area)
              if gamma(ii+1)>=gamma(ii-1)
                  Aii = [A_theta(:,ii) A_theta(:,ii+1)];
                  gammaii = diag([gamma(ii);gamma(ii+1)]);
              else
                  Aii = [A_theta(:,ii-1) A_theta(:,ii)];
                  gammaii = diag([gamma(ii-1);gamma(ii)]);                 
              end
              sigmaYii = sigmaY - Aii*gammaii*Aii';
              
              %calculate zm qm
              search_gap = 0.1;
              search_temp = search_area(ii-1):search_gap:search_area(ii+1);
              A_temp = exp(-1i*2*pi*(0:M-1)'*(dd/lambda)*cosd(search_temp));
              Zm = diag(A_temp'/sigmaYii*A_temp);
              qm_temp = abs(A_temp'/sigmaYii*Y).^2;
              qm = sum(qm_temp,2);
              
              %calculate gammam
              e = (scale_lambda - 1 - T).*Zm.^2;
              f = (2*scale_lambda - 2 - T).*Zm + qm;
              g = scale_lambda - 1;
              gammam = (-f - sqrt(f.^2 - 4*g*e))./(2*e);
              
              %maximize L
              L = -T*log(1+gammam.*Zm) + qm./(1./gammam+Zm) + g*log(gammam);
              [~,indexmax] = max(L);
              search_area(ii) =  search_temp(indexmax);                         
          end
      end
      
Pm=mean(mu.*conj(mu),2);

end

