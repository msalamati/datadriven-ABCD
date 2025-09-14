clc
clear all
close all
%%
% Seting the example-specific hyper-parameters
num_points = 20; % number of points for plotting
beta = .01; % confidence = 1-beta
eps =.01; % this constant affects the level of conservativeness
eta_x = .2; % state discretization size
eta_u = .3; % input discretization size
w_bar = .01; % upper-bound on the value of disturbance vector
dim_x = 3; % dimension of the state space
dim_u = 2; % dimension of the state space
L_d = 1.46; % upper-bound over the Lipschitz constant 
is_adaptive = 0; % set this binary variable to True for using the adaptive GB
dist_presence = 1; % set this binary variable to True for sampling the disturbance space
if dist_presence
    n = dim_x*2; % sample both the state and disturbance spaces
else
    n = dim_x; % only sample the state space
end
if is_adaptive
    n = n*2; % sample both the state and disturbance spaces TWICE
end
N_x = ceil(10/eta_x)*ceil(10/eta_x)*ceil((2*pi+.8)/eta_x); % number of discrete states (for examples other than unicycle pls change it)
N_u = ceil(1/eta_u)*ceil(1/eta_u); % number of discrete states (for examples other than unicycle pls change it)
conf_per_cell = beta;%/N_x*N_u; % modify the overall confidence (1-conf_per_cell) based on number of cells
%%
% plot effect of changing epsilon on sample complexity
figure()
N=zeros(num_points,1);
epss =zeros(num_points,1);
for i=1:1:num_points
    epss(i)= eps*i;
    N(i) = n-1; % initializing sample numbers
    lhs =lhs_computer(N(i),n,epss(i));
    while lhs>conf_per_cell
        N(i) =N(i)+1;
        lhs =lhs_computer(N(i),n,epss(i));
    end
    if i ==1
        N_max=N(i);
    else
        N_max=max(N(i),N_max);
    end 
end
plot(epss, N,'-*','color', 'r', 'LineWidth',2)
hold on
N1=max(N);
N=zeros(num_points,1);
epss =zeros(num_points,1);
for i=1:1:num_points
    epss(i)= eps*i;
    N(i) = dim_x-1; % initializing sample numbers
    lhs =lhs_computer(N(i),dim_x,epss(i));
    while lhs>conf_per_cell
        N(i) =N(i)+1;
        lhs =lhs_computer(N(i),dim_x,epss(i));
    end
    if i ==1
        N_max=N(i);
    else
        N_max=max(N(i),N_max);
    end 
end
plot(epss, N,'-*','color', 'b', 'LineWidth',2)
N2=max(N);
legend('With disturbance','Without disturbance', 'interpreter','latex', 'FontSize', 20)
%legend('Variations of $N$ with $\varepsilon$ when $W\neq\{0\}$','Sample number in the absence of disturbance', 'interpreter','latex', 'FontSize', 15)

xlabel('$\varepsilon$','interpreter','latex', 'FontSize', 20)
ylabel('$N$','interpreter','latex', 'FontSize', 20)
xlim([eps,1.01*eps*num_points])
ylim([0,1.01*max(N1,N2)])
%title('Effect of Variations in required ')
grid
set(gca,'FontSize',20)
axis tight
%%
% plot effect of changing beta on sample complexity
N=zeros(num_points,1);
conf_per_celll =zeros(num_points,1);
figure()
for i=1:1:num_points
    conf_per_celll(i) = beta*i;
    N(i) = n-1; % initializing sample numbers
    lhs =lhs_computer(N(i),n,eps);
    while lhs>conf_per_celll(i)
        N(i) =N(i)+1;
        lhs =lhs_computer(N(i),n,eps);
    end
    if i ==1
        N_max=N(i);
    else
        N_max=max(N(i),N_max);
    end
end
plot(conf_per_celll, N,'-*','color', 'r', 'LineWidth',2)
hold on
N1=max(N);
for i=1:1:num_points
    conf_per_celll(i) = beta*i;
    N(i) = dim_x-1; % initializing sample numbers
    lhs =lhs_computer(N(i),dim_x,eps);
    while lhs>conf_per_celll(i)
        N(i) =N(i)+1;
        lhs =lhs_computer(N(i),dim_x,eps);
    end
    if i ==1
        N_max=N(i);
    else
        N_max=max(N(i),N_max);
    end
end
plot(conf_per_celll, N,'-*','color', 'b', 'LineWidth',2)
N2=max(N);
legend('With disturbance','Without disturbance', 'FontSize', 20)
xlabel('$\beta$','interpreter','latex', 'FontSize', 20)
ylabel('$N$','interpreter','latex', 'FontSize', 20)
xlim([beta,1.01*beta*num_points])
ylim([0,1.01*max(N1,N2)])
grid
set(gca,'FontSize',20)
axis tight
%%
% plot effect of changing epsilon on bias value
bias=zeros(num_points,1);
epss =zeros(num_points,1);
figure()
for i=1:1:num_points
    epss(i)= eps*i;
    %L_d = max(eta_x/2,1); % the value of Lipschitz constant in theta
    bias(i) = L_d*(((eta_x/2)^dim_x)*((2*w_bar)^dim_x)*epss(i))^(1/n); % the bias
end
plot(epss, bias, '-*', 'color', 'r', 'LineWidth',2)
hold on
b1=max(bias);

for i=1:1:num_points
    epss(i)= eps*i;
    %L_d = max(eta_x/2,1); % the value of Lipschitz constant in theta
    bias(i) = L_d*(((eta_x/2)^dim_x)*epss(i))^(1/dim_x); % the bias
end
plot(epss, bias, '-*', 'color', 'b', 'LineWidth',2)
b2=max(bias);
legend('With disturbance','Without disturbance', 'FontSize', 20)
xlabel('$\varepsilon$','interpreter','latex', 'FontSize', 20)
ylabel('$\gamma$','interpreter','latex','FontSize', 20)
xlim([eps,1.01*eps*num_points])
ylim([0,1.01*max(b1,b2)])
grid
set(gca,'FontSize',20)
axis tight
%%
% plot effect of changing eta_x on bias value...ETA_X HAS NO EFFECT ON THE
% BIAS VALUE
% bias=zeros(num_points,1);
% eta_xx =zeros(num_points,1);
% figure()
% for i=1:1:num_points
%     eta_xx(i)= eta_x*i;
%     L_d = max(eta_xx(i)/2,1); % the value of Lipschitz constant in theta
%     bias(i) = L_d*((((eta_xx(i)/2)^dim_x)*((2*w_bar)^dim_x)*eps)^(1/n)); % the bias
% end
% plot(eta_xx, bias, '-*', 'color', 'r', 'LineWidth',2)
% hold on
% b1 = max(bias);
% for i=1:1:num_points
%     eta_xx(i)= eta_x*i;
%     L_d = max(eta_xx(i)/2,1); % the value of Lipschitz constant in theta
%     bias(i) = L_d*(((eta_xx(i)/2)^dim_x)^(1/dim_x)); % the bias
% end
% plot(eta_xx, bias, '-*', 'color', 'b', 'LineWidth',2)
% hold on
% b2 = max(bias);
% legend('With disturbance','Without disturbance', 'FontSize', 20)
% xlabel('$\eta_x$','interpreter','latex', 'FontSize', 20)
% ylabel('$\gamma$','interpreter','latex','FontSize', 20)
% xlim([eta_x,1.01*eta_x*num_points])
% ylim([0,1.01*max(b1,b2)])
% grid
% 
% set(gca,'FontSize',20)
% axis tight


%%
% defining the necessary functions
function lhs = lhs_computer(N,n,eps)
    lhs = 0;
    for i=0:1:n-1
        lhs = lhs+nchoosek(N,i)*eps^i*(1-eps)^(N-i);
    end
end


