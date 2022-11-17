function [Y] = func_AR1_2(phi, T, B, sigma)
nb_process = B;

Y=nan(T,nb_process);
Y(1,:)=zeros(1,nb_process); % First row is zero

for j=1 : nb_process
    for t=2 : T
        Y(t,j) = phi * Y(t-1, j) + randn() * sigma;
    end
end
