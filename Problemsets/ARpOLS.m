%%%% OLS Function %%%%

%%% Inputs %%%
% y : data vectors
% p : number of lags
% const = 1 if constant, or 2 if constant + linear trend
% alpha : significance of the p-values

%function [thetas, std_errs, t_stats, p_values] = ARpOLS(Y,p,const,alpha)
function [thetas, SE_thetas, t_stats, y, y_hat] = ARpOLS(Y,p,const)
len = length(Y)-p;
time_vect = nan(len,1);
const_vect = zeros(len, 1)+1;
y = Y(1+p:end, 1);
for i = 1 : len 
    time_vect(i,1)=i;
end

if const ==1 
    y_lagged = [const_vect Y(1:end-p, 1)];
elseif const ==2
    y_lagged = [const_vect time_vect Y(1:end-p, 1)];
end

thetas = ((y_lagged'*y_lagged)^(-1))*(y_lagged'*y);

y_hat = sum((Y(1:end-p, 1)*thetas'), 2);
u_hat =  y-y_hat;

sigma_square_u=((1/(len-1-p))*sum(u_hat));
var_thetas=sigma_square_u*((y_lagged'*y_lagged)^(-1));
SE_thetas = diag(sqrt(var_thetas));

t_stats = (thetas-1)/SE_thetas;


end