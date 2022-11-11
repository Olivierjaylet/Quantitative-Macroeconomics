%%%% OLS Function %%%%

%%% Inputs %%%
% y : [Tx1] data vectors
% p : [scalar] number of lags
% const = 1 if constant, or 2 if constant + linear trend
% alpha : significance of the p-values

function OLS = ARpOLS(Y,p,const, alph)
len = size(Y,1)-p;
time_vect = nan(len,1);
const_vect = zeros(len, 1)+1;
y = Y(1+p:end, 1);
for i = 1 : len 
    time_vect(i,1)=i;
end

if const ==1 
    y_lagged = [const_vect Y(1:end-p, 1)];
elseif const == 2
    y_lagged = [const_vect time_vect Y(1:end-p, 1)];
end

YtYinv = inv(y_lagged'*y_lagged);
thetahat = YtYinv*(y_lagged'*y);

yhat = sum((Y(1:end-p, 1)*thetahat'), 2);
uhat =  y-yhat;
utu = uhat'*uhat;

var_uhat=utu/(len-p-const);
siguhat = sqrt(var_uhat);

var_thetahat=var_uhat*(diag(YtYinv));
sd_thetahat = sqrt(var_thetahat);

t_stats = thetahat./sd_thetahat;
t_crit = -tinv(alph/2,len-p-const);
pvalues = tpdf(t_stats,len-p-const);

theta_ci = [thetahat-t_crit.*sd_thetahat, thetahat+t_crit.*sd_thetahat];

OLS.T_eff = len;
OLS.thetahat = thetahat;
OLS.siguhat = siguhat;
OLS.sd_thetahat = sd_thetahat;
OLS.tstats = t_stats;
OLS.pvalues = pvalues;
OLS.theta_ci = theta_ci;
OLS.resid = uhat;

end