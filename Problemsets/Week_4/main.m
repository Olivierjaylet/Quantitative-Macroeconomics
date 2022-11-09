Y=func_AR1_2(0.95, 1000, 1, 1);

[thetas, SE_thetas, t_stats, y, y_hat] = ARpOLS(Y,1,1);
figure;
scatter(y_hat, y, 'blue');
a= lsline;
set(a(1),'color','r');
