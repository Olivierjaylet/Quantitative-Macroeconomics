%%%% Conditional log-likelihood contributions %%%%

function loglik = func_CondLogLike(x, Y, p, const)

T=size(Y);
y = Y(1+p:end, 1); % Y is the lagged matrix and y at t-time

theta = x(1:(const+p)); % All values before sig_u
sig_u = x(const+p+1); % Last value of x

if const == 1
    Y = [ones(T,1) Y];
elseif const ==2
    Y = [ones(T,1) transpose(1:T) Y];
end

% Compute residuals
uhat = y - Y*theta;
% Compute Sum of squared residual
SSR = transpose(uhat)*uhat;

loglik = -log(2*pi)*(T-p)/2 - log(sig_u^2)*(T-p)/2 - SSR/(2*sig_u^2);


% c, fi and S are in the parameters vector
c = parameters(1,:); % first row of the vector
fi = parameters(2,:); % second
s = parameters(3,:); % third

Yt = timeseries(2:end, 1); % value of y at t time
Y_1 = timeseries(1: end-1, 1); % value of y at t-1 time
eps = Yt - c - fi * Y_1; % calculate epsilon, the sequence vector 

condloglike=[];
for i=2:T-1
    condloglike(i,:) = log(1/sqrt(2*pi*(s)))-((eps(i,:)^2)/(2*(s)));
end

end
