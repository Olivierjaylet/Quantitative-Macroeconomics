%%%% Conditional log-likelihood contributions %%%%

function loglik = func_LogLikeARpNorm(x, y, p, const)

T=size(y);
Y = lagmatrix(y,1:p);

theta = x(1:(const+p)); % AR coefficients (All values before sig_u)
sig_u = x(const+p+1);   % Last value of x

if const == 1
    Y = [ones(T,1) Y];
elseif const ==2
    Y = [ones(T,1) transpose(1:T) Y];
end
Y = Y((p+1):end,:);     % Get rid of initial observations
y = y(p+1:end);         % Get rid of initial observations

% Compute residuals
uhat = y - Y*theta;
% Compute Sum of squared residual
SSR = transpose(uhat)*uhat;


% Compute the conditional log likelihood
loglik = -log(2*pi)*(T-p)/2 - log(sig_u^2)*(T-p)/2 - SSR/(2*sig_u^2);

if isnan(loglik) || isinf(loglik) || ~isreal(loglik)
    loglik = -1e10;     % if anything goes wrong set value very small, can also use -Inf
end

end