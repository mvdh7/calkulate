    

TD = [288.15  1.029691
      290.65  1.029084
      293.15  1.028426
      295.65  1.027720
      298.15  1.026966
      300.65  1.026169
      303.15  1.025328
      305.65  1.024446
      308.15  1.023524];
 
T = TD(:,1);
D = TD(:,2);

load('HCl in NaCl.mat')

fx = 278.15:0.05:318.15;
fy = polyval(pfit,fx);

pfit2 = [-3.59047619e-06,  1.83214095e-03,  7.99883606e-01];

fy2 = polyval(pfit2,fx);

fy3 = pfit2(1)*fx.^2 + pfit2(2)*fx + pfit2(3);

figure(1); clf; hold on

scatter(TD(:,1),TD(:,2))

plot(fx,fy)
plot(fx,fy2)
plot(fx,fy3)
