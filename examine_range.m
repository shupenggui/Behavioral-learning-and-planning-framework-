function [range] = examine_range(data,Y,lambda,iteration)
	a = size(Y,2);
    b = size(data,2);
	W = zeros(a,b);
    W(1,:) = 1;
	V = W;
	u = zeros(a,b);
	rho =2;
	drho = 1/rho;
	beta = lambda*drho;
%     obj = zeros(1,iteration);
	D = Y'*data;
	Dt = Y * (Y');
    delta = inv(eye(4)+ drho * Dt);
    Ai1 = drho * ones(a,1) - (drho^2 * (Y')) * ((delta) * (Y * ones(a,1)));
    N1 = Ai1/(ones(1,a)*Ai1);
    
    for i = 1:iteration
%         obj(i) = 0.5*norm(data-Y*W,'fro')^2 + lambda * sum(sqrt(sum(W.^2,2)));
	    B = D+rho*V - u;
        W = drho * B - (drho^2 * (Y')) * (((delta) * Y )*B) - N1*(Ai1'*B - ones(1,b));
        M = max(0,W + drho * u);
	    V = max(0,1-beta./(sqrt(sum(M.^2,2))+eps))*ones(1,b).*M;
	    u = u + rho*(W-V);
    end
    
% 	figure;
% 	semilogy(obj);
    range = W;
    

end