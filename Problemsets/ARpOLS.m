%%%% OLS Function %%%%

%%% Input %%%
% y : [Tx1] data vectors
% p : [scalar] number of lags
% const = 1 if constant, or 2 if constant + linear trend
% alpha : significance of the p-values

%%% Output %%%
% OLS Structure :
%   - T_eff :       [scalar]        effective sample size used in estimation
%   - thetathat :   [(const+p)x1]   estimation of coefficients
%   - siguhat :     [scalar]        estimate of standard deviation of error term
%   - sd_thetahat : [(const+p)x1]   estimate of standard error of coefficients
%   - tstats :      [(const+p)x1]   t statistics
%   - pvalues :     [(const+p)x1]   p values of H_0 : thetahat=0
%   - theta_ci :    [(const+p)x2]   (1-alph)% confidence intervall for theta given
%                                   significance level alph
%   - resid :       [T_effx1]       residuals



function OLS = ARpOLS(y,p,const, alph)
T=size(y,1);                            % sample size
T_eff = T-p;                            % effective sample size used in estimation
time_vect = transpose(1:T);     % time term
const_vect = ones(T, 1);            % constant term
Y = lagmatrix(y,1:p);                   % create matrix xith lagged variables


if const ==1
    Y = [const_vect Y];             % constant term
elseif const == 2
    Y= [const_vect time_vect Y];    % constant term + time trend
end

Y=Y((p+1):end,:); % get rid of initial p observations
y=y((p+1):end,:); % get rid of initial p observations

YtYinv = inv(Y'*Y);
thetahat = YtYinv*(Y'*y);   % OLS estimator of coefficients;

yhat = Y*thetahat;          % predicted values
uhat =  y-yhat;             % residuals
utu = uhat'*uhat;           

var_uhat=utu/(T_eff-p-const);           % variance of error term
siguhat = sqrt(var_uhat);               % standard deviation of error term
var_thetahat=var_uhat*(diag(YtYinv));   % variance of coefficients
sd_thetahat = sqrt(var_thetahat);       % sd of coefficients

t_stats = thetahat./sd_thetahat;        % t-statistics
t_crit = -tinv(alph/2,T_eff-p-const);     % critical value
pvalues = tpdf(t_stats,T_eff-p-const);    % p-values

% confidence interval
theta_ci = [thetahat-t_crit.*sd_thetahat, thetahat+t_crit.*sd_thetahat];



OLS.T_eff = T_eff;
OLS.thetahat = thetahat;
OLS.siguhat = siguhat;
OLS.sd_thetahat = sd_thetahat;
OLS.tstats = t_stats;
OLS.pvalues = pvalues;
OLS.theta_ci = theta_ci;
OLS.resid = uhat;

end