function [Pm,iter,Error,Gamma_a] = SBL_gh(Y,dd,lambda,search_area)

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
Error = zeros(maxiter,1);
Gamma_a = zeros(N,maxiter);

while ~converged
      
      iter = iter + 1;
      gamma_last = gamma;
      
      %calculate mu and Sigma
      sigmaY = eye(M)/beta + A_theta*diag(gamma)*A_theta';
      mu = diag(gamma)*A_theta'/sigmaY*Y;
      sigma = diag(gamma) - diag(gamma)*A_theta'/sigmaY*A_theta*diag(gamma);
      
      %calculate gamma
      gamma_b = gamma_b0 + 2*diag(mu*mu'+ T*sigma);
      gamma_lambda = gamma_lambda0 - T;
      gamma_a = (k_a1+gamma_lambda0/2)./(k_a2 + gamma/2);
      gamma = sqrt(gamma_b./gamma_a).* besselk(gamma_lambda+1,sqrt(gamma_a.*gamma_b))./besselk(gamma_lambda,sqrt(gamma_a.*gamma_b));
      Gamma_a(:,iter) = gamma_a;
      
      %update beta
      resid = Y-A_theta*mu;
      beta = (M*T+c)/(norm(resid, 'fro')^2+T*trace(A_theta*sigma*A_theta')+d);

      %stopping criteria
      erro=norm(gamma - gamma_last)/norm(gamma_last);
      Error(iter) = erro;
      if erro < tol || iter >= maxiter
          converged = true;
      end
      
end
      
Pm=mean(mu.*conj(mu),2);

end

