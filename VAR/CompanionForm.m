function A = CompanionForm(Coefs, p)

% Computes matrix A of companion form of a VAR(p) model with constant

K = size(Coefs, 1); % number of variables
const = size(Coefs, 2) - K*p;

A = [Coefs(:,const+1:end);
    eye(K*(p-1)) zeros(K*(p-1), K)];

end