# This file will generate the plot about ARE between T_{tW} and T_{W}. 
# We consider 4 cases here. 
# The result figure is called "ARE_Wilcoxon.pdf"

pdf(file = "ARE_Wilcoxon.pdf", width = 8, height = 8)
par(mfrow = c(2, 2))

N1 = 40; N2 = 50;
theta1 = 0.3; theta2 = 0.8;
n1 = N1*theta1; n2 = N2*theta2;
thetam = max(theta1, theta2);
theta = (theta1 + theta2)/2;
dtheta = theta2 - theta1;
d12 = seq(-0.5, 0.5, 0.01);

vS = N1*N2*(N1+N2)*theta*(theta^2+3*(1 - theta))/12;
vs = N1*N2*(N1+N2)*theta^3*(1+3*(1 - theta))/12;
meanS = n1*n2*d12 + N1*N2*dtheta/2;
means = n1*n2*d12 + N1*N2*dtheta/2*thetam - theta*(1 - theta)*(N2 - N1)/4;
ratio = (means/meanS)^2*(vS/vs);

par(mar=c(4,4,2,2)) 
plot(d12, ratio, type = "l", lwd = 2, xlab = expression(Delta["f,g"]), 
ylab = expression("ARE(T"[tW]*", T"[W]*")"))
text(x = -0.48, y = 1.65, "(a)", pos = 3, cex = 1.2)
abline(h = 1, lty = 2)
# text(-0.48, 1.05, "y=1", cex = 1)

N1 = 40; N2 = 50;
theta = seq(0.1, 0.9, by = 0.01);
theta1 = theta - 0.1; theta2 = theta + 0.1;
n1 = N1*theta1; n2 = N2*theta2;
thetam = apply(cbind(theta1, theta2), 1, max);
theta = (theta1 + theta2)/2;
dtheta = theta2 - theta1;
d12 = 0.1;

vS = N1*N2*(N1+N2)*theta*(theta^2+3*(1 - theta))/12;
vs = N1*N2*(N1+N2)*theta^3*(1+3*(1 - theta))/12;
meanS = n1*n2*d12 + N1*N2*dtheta/2;
means = n1*n2*d12 + N1*N2*dtheta/2*thetam - theta*(1 - theta)*(N2 - N1)/4;
ratio = (means/meanS)^2*(vS/vs);

par(mar=c(4,4,2,4)) 
plot(theta, ratio, type = "l", lwd = 2, xlab = expression(theta), ylab = expression("ARE(T"[tW]*", T"[W]*")"), 
     ylim = c(1, max(ratio)))
text(0.85, 2.65, "(b)", pos = 3, cex = 1.2)
abline(h = 1, lty = 2)
# text(0.15, 1.05, "y=1", cex = 1)

N1 = 40; N2 = 50;
dtheta = seq(0, 0.5, by = 0.01);
theta1 = 0.5 - dtheta; theta2 = 0.5 + dtheta;
n1 = N1*theta1; n2 = N2*theta2;
thetam = apply(cbind(theta1, theta2), 1, max);
theta = (theta1 + theta2)/2;
d12 = 0.1;

vS = N1*N2*(N1+N2)*theta*(theta^2+3*(1 - theta))/12;
vs = N1*N2*(N1+N2)*theta^3*(1+3*(1 - theta))/12;
meanS = n1*n2*d12 + N1*N2*dtheta/2;
means = n1*n2*d12 + N1*N2*dtheta/2*thetam - theta*(1 - theta)*(N2 - N1)/4;
ratio = (means/meanS)^2*(vS/vs);

par(mar=c(4,4,2,2)) 
plot(dtheta*2, ratio, type = "l", lwd = 2, xlab = expression(Delta*theta), ylab = expression("ARE(T"[tW]*", T"[W]*")"), 
     ylim = c(1, max(ratio)))
text(0.03, 2.65, "(c)", pos = 3, cex = 1.2)
abline(h = 1, lty = 2)
# text(0.1, 1.05, "y=1", cex = 1)

N1 = seq(20, 120, by = 5); N2 = 50;
theta1 = 0.3; theta2 = 0.8;
n1 = N1*theta1; n2 = N2*theta2;
thetam = max(theta1, theta2);
theta = (theta1 + theta2)/2;
dtheta = theta1 - theta2;
d12 = 0.1;

vS = N1*N2*(N1+N2)*theta*(theta^2+3*(1 - theta))/12;
vs = N1*N2*(N1+N2)*theta^3*(1+3*(1 - theta))/12;
meanS = n1*n2*d12 + N1*N2*dtheta/2;
means = n1*n2*d12 + N1*N2*dtheta/2*thetam - theta*(1 - theta)*(N2 - N1)/4;
ratio = (means/meanS)^2*(vS/vs);

ratio = (1/theta^2)*(1 - theta + theta^2/3)/(1 - theta + 1/3)*((theta1*theta2*d12 + dtheta*thetam/2)/(theta1*theta2*d12 + dtheta/2))^2
par(mar=c(4,4,2,4)) 
plot(N1, rep(ratio, length(N1)), type = "l", lwd = 2, xlab = expression('N'[1]), ylab = expression("ARE(T"[tW]*", T"[W]*")"))
text(116, 1.8, "(d)", pos = 3, cex = 1.2)
abline(h = 1, lty = 2)
# text(min(N1) + 2, 1.05, "y=1", cex = 1)
dev.off()

