function [Y] = func_AR1_1(phi, T, sigma)
nb_process = length(phi);

Y=nan(T,nb_process);
Y(1,:)=zeros(1,nb_process); % First row is zero

for j=1 : nb_process
    for t=2 : T
        Y(t,j) = phi(j) * Y(t-1, j) + randn() * sigma;
    end
end
