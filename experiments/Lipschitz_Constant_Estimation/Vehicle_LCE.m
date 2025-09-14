%author: B.Wooding
%date: Jan 2022
clc
clear all
close all
%% Algorithm
% Calculating an estimate of the Lipschitz Constant M
%% Step 1: Sample the slopes
% Given delta > 0, choose pairs (x_i,y_i) uniformly on
%{(x,y) \elems [a,b] x [a,b] : |x-y| <= delta}
% and evaluate
% s_i = |g(x_i) - g(y_i)|/|x_i-y_i| for i = 1,...,n

%% Step 2: Calculate the Maximum Slope
% l = max{s_i,...,s_n}

% Steps 1 and 2 are performed m times, giving l_1, ..., l_m.

%% Step 3: Fit the Reverse Weibull Distribution (RWD)
% Fit a three parameter Reverse Weibull distribution to l_1, ..., l_m.

%% Output:

% Our estimate Mhat, of M, is the location parameter of the fitted Reverse
% Weibull distribution.

% It is remarked that the Lipschitz constant estimation problem is converted into a
% routine curve fitting problem.

%% Source
% Wood, G.R., Zhang, B.P. Estimation of the Lipschitz constant of a function.
% J Glob Optim 8, 91–103 (1996). https://doi.org/10.1007/BF00229304

%% Further Remarks
% For our system additional external loops need to be added to have a
% Lipschitz Constant estimated for each input and each disturbance. We take
% the discrete inputs and disturbances used in the SCOTS abstraction for
% these.

%% Parameters
n = 100; % number of slopes
m = 1000; % number of data points for RWD
dim = 3; % state dimensions
% list of all inputs
u1 = [-1; -0.9; -0.6; -0.3; -0; 0.3; 0.6; 0.9; 1];
u2 = [-1; -0.9; -0.6; -0.3; -0; 0.3; 0.6; 0.9; 1];
i_dim = 2; % input dimensions
ins1 = size(u1,1); % number of inputs in dim 1
ins2 = size(u2,1); % number of inputs in dim 2
M = zeros(ins1*ins2,1);  % store lipschitz constant estimations for each input
lb = [0  0 -pi-0.4]; %lowerbound of system domain
ub = [10 10 pi+0.4]; %upperbound of system domain

% uniform distributions [a,b] x [a,b] for x_i and y_i
dist = {
    makedist('Uniform','lower',lb(1),'upper',ub(1));
    makedist('Uniform','lower',lb(2),'upper',ub(2));
    makedist('Uniform','lower',lb(3),'upper',ub(3));
    };
%% Implementation
% loop through all inputs
for iy = 1:ins2
    for ix = 1:ins1
        %data store for RWD
        l = zeros(m,1);
        % loop through to calculate data for RWD
        for k = 1:m
            % x_i and y_i pairs for slope
            xx = zeros(dim,1);
            yy = zeros(dim,1);
            slopes = zeros(n,1);
            %% Step 1: Sample the slopes
            % loop through slope calculations
            for j = 1:n
                for i=1:dim
                    % Choose random x_i and y_i
                    xx(i,:)=random(dist{i});
                    yy(i,:)=random(dist{i});
                end
                %calculate s_i
                u = [u1(ix) u2(iy)];
                [t1 ode1]=ode45(@g,[0 .4], xx, [],u);
            [t2 ode2]=ode45(@g,[0 .4], yy, [],u);
            %g(xx,u(ix))
            %g(yy,u(ix))
            slopes(j) = norm(ode1(end,:) - ode2(end,:))/norm(xx-yy);
            end
            %% Step 2: Calculate the Maximum Slope
            %Store l as data for RWD
            l(k) = max(slopes);
        end
        %% Step 3: Fit a three-parameter Reverse Weibull distribution
        %l = data
        %v = scale
        %w = shape
        %u = location parameter
        
        %custom pdf for the three parameter reverse weibull
        f = @(l,v,w,u) (l<=u).*(w.*(u-l).^(w-1).*exp((-(u-l).^w)/v)/v);
        
        %fit data and calculate parameters
        opt = statset('FunValCheck','off');
        params = mle(l,'pdf',f,'Start',[0.76482 5.61034 max(l)],'Options',opt,'LowerBound',[0 0 max(l)],'UpperBound',[Inf Inf Inf]);
        M(ins2*(iy-1) + ix) = params(3); % save location parameter to set of Lipschitz Constants
    end
    iy
end
%% Output
max(M)
%% functions
function dxdt = g(t,x,u)
alpha=atan(tan(u(2))/2.0);
dxdt(1) = u(1)*cos(alpha+x(3))/cos(alpha);
dxdt(2) = u(1)*sin(alpha+x(3))/cos(alpha);
dxdt(3) = u(1)*tan(u(2));

dxdt = [dxdt(1); dxdt(2); dxdt(3)];
end