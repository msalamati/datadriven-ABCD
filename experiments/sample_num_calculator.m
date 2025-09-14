clc
clear all
close all
%%
% Seting the example-specific hyper-parameters
beta = .01; % confidence = 1-beta
eps =.01; % this constant affects the level of conservativeness
is_adaptive = 0; % set this binary variable to True for using the adaptive GB
dist_presence = 0; % set this binary variable to True for sampling the disturbance space

%3d vehicle
L_d = 1.3; % the value of Lipschitz constant in uncertainty
eta_x = .2; % state discretization size
eta_u = .3; % input discretization size
w_bar = .005; % upper-bound on the value of disturbance vector
dim_x = 3; % dimension of the state space
dim_u = 2; % dimension of the state space
N_x = (ceil(5/eta_x)+1)*(ceil(5/eta_x)+1)*(ceil((2.2)+1)/eta_x); % number of discrete states (for examples other than unicycle pls change it)
N_u = ceil(2.2/eta_u)*ceil(2.2/eta_u); % number of discrete states (for examples other than unicycle pls change it)

%4d vehicle
% L_d = 2; % the value of Lipschitz constant in uncertainty
% eta_x = .2; % state discretization size
% eta_u = .3; % input discretization size
% w_bar = .005; % upper-bound on the value of disturbance vector
% dim_x = 4; % dimension of the state space
% dim_u = 2; % dimension of the state space
% N_x = (ceil(5/eta_x)+1)*(ceil(5/eta_x)+1)*(ceil((2.2)+1)/eta_x)*(ceil((2.2)+1)/eta_x); % number of discrete states (for examples other than unicycle pls change it)
% N_u = ceil(2.2/eta_u)*ceil(2.2/eta_u); % number of discrete states (for examples other than unicycle pls change it)

%5d vehicle
% L_d = 2; % the value of Lipschitz constant in uncertainty
% eta_x = .2; % state discretization size
% eta_u = .3; % input discretization size
% w_bar = .005; % upper-bound on the value of disturbance vector
% dim_x = 5; % dimension of the state space
% dim_u = 2; % dimension of the state space
% N_x = (ceil(5/eta_x)+1)*(ceil(5/eta_x)+1)*(ceil((2.2)+1)/eta_x)*(ceil((2.2)+1)/eta_x)*(ceil((2.2)+1)/eta_x); % number of discrete states (for examples other than unicycle pls change it)
% N_u = ceil(2.2/eta_u)*ceil(2.2/eta_u); % number of discrete states (for examples other than unicycle pls change it)


%2d dcdc
% L_d = 1.3; % the value of Lipschitz constant in uncertainty
% eta_x = 20/4000; % state discretization size
% eta_u = 1; % input discretization size
% w_bar = .005; % upper-bound on the value of disturbance vector
% dim_x = 2; % dimension of the state space
% dim_u = 1; % dimension of the state space
% N_x = (ceil(1/eta_x)+1)*(ceil(1/eta_x)+1); % number of discrete states 
% N_u = ceil(1/eta_u)+1; % number of discrete states


if dist_presence
    n = dim_x*2; % sample both the state and disturbance spaces
else
    n = dim_x; % only sample the state space
end
if is_adaptive
    n = n*2; % sample both the state and disturbance spaces TWICE
end
conf_per_cell = beta/(N_x*N_u); % modify the overall confidence (1-conf_per_cell) based on number of cells
%%
% find the minimum number of samples
N = n-1; % initializing sample numbers
lhs =lhs_computer(N,n,eps);
while lhs>conf_per_cell
    N =N+1;
    lhs =lhs_computer(N,n,eps);
end
to_be_displayed1 = ['minimum number to be sampled for each cell is ', num2str(N)];
disp(to_be_displayed1);
%%
% compute the value of the bias that must be added to the computed GB
% L_d = max(eta_x/2,w_bar); % the value of Lipschitz constant in theta
if is_adaptive
    if dist_presence
        bias = L_d*(((eta_x/2)^(2*dim_x))*((2*w_bar)^(2*dim_x))*eps)^(1/n);
    else
        bias = L_d*(((eta_x/2)^(2*dim_x))*eps)^(1/n);
    end
else
    if dist_presence
        bias = L_d*(((eta_x/2)^dim_x)*((2*w_bar)^dim_x)*eps)^(1/n);
    else
        bias = L_d*(((eta_x/2)^dim_x)*eps)^(1/n);
    end
end
to_be_displayed2 = ['The bias which should be added to the data-driven GB is ', num2str(bias)];
disp(to_be_displayed2);


%%
% computing the sample number using the methot in Murat's paper
eps = 0.1;
N = ceil((N_x*log(2)+log((1-eps)/(N_x*N_u)))/(1-beta));
to_be_displayed3 = ['Sample number using the method in Murat paper is ', num2str(N)];
disp(to_be_displayed3);


%%
% defining the necessary functions
function lhs = lhs_computer(N,n,eps)
    lhs = 0;
    for i=0:1:n-1
        lhs = lhs+nchoosek(N,i)*eps^i*(1-eps)^(N-i);
    end
end


