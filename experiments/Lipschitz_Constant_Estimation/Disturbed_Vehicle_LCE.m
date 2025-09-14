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
rng(1); % set random seed to 1
%% Parameters
n = 3; % number of slopes
m = 5; % number of data points for RWD
dim = 3; % state dimensions
% list of all inputs
u1 = [-0.6; -0.3; -0; 0.3; 0.6]; %include 0.6,0.9, 1, -1, -0.9 -0.6
u2 = [-0.6; -0.3; -0; 0.3; 0.6]; %include 0.6,0.9, 1, -1, -0.9 -0.6
i_dim = 2; % input dimensions
ins1 = size(u1,1); % number of inputs in dim 1
ins2 = size(u2,1); % number of inputs in dim 2
M = zeros(ins1*ins2,1);  % store lipschitz constant estimations for each input
eta = 0.5; %10*2/4e3;
dim1 = [0:eta:10];
dim2 = [0:eta:10];
[p,q] = meshgrid(dim1, dim2);
pairs = [p(:) q(:)];
LSol = zeros(size(pairs,1),1);
%lb = [0  0 -pi-0.4]; %lowerbound of system domain
%ub = [10 10 pi+0.4]; %upperbound of system domain

for pair = 1:size(pairs,1)
    % uniform distributions [a,b] x [a,b] for x_i and y_i
    dist = {
        makedist('Uniform','lower',pairs(pair,1),'upper',pairs(pair,1)+eta);
        makedist('Uniform','lower',pairs(pair,2),'upper',pairs(pair,2)+eta);
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
                    for i=1:dim-1
                        % Choose random x_i and y_i
                        xx(i,:)=random(dist{i});
                        yy(i,:)=random(dist{i});
                    end
                    %calculate s_i
                    u = [u1(ix) u2(iy)];
                    [t1 ode1]=ode45(@g,[0 .5], xx, [],u);
                    [t2 ode2]=ode45(@g,[0 .5], yy, [],u);
                    %g(xx,u)
                    %g(yy,u)
                    slopes(j) = norm(ode1(end,:) - ode2(end,:))/norm(xx-yy);
                    %slopes(j) = norm(g(xx,u) - g(yy,u))/norm(xx-yy);
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
            %pd = fitdist(l,'Weibull')
            %fit data and calculate parameters
            opt = statset('FunValCheck','off');
            % 0.76482 5.61034
            params = mle(l,'pdf',f,'Start',[0.52 9.69 max(l)],'Options',opt,'LowerBound',[0 0 max(l)],'UpperBound',[Inf Inf Inf]);
            %0.52 9.69 max(l)
            M(ins2*(iy-1) + ix) = params(3); % save location parameter to set of Lipschitz Constants
        end
    end
    %% Output
    LSol(pair) = max(M);
    pair
end
%%
%[xq,yq] = meshgrid(-10:.2:10, -10:.2:10);
vq = griddata(pairs(:,1),pairs(:,2), LSol,p,q);  %(x,y,v) being your original data for plotting points
mesh(p,q,vq)
hold on
plot3(pairs(:,1),pairs(:,2), LSol,'o')
zlim([1.5 1.6])

%% functions
function dxdt = g(t,x,u)
w = 0.01;
alpha=atan(tan(u(2))/2.0);
dxdt(1) = u(1)*cos(alpha+x(3))/cos(alpha) + w;
dxdt(2) = u(1)*sin(alpha+x(3))/cos(alpha);
dxdt(3) = u(1)*tan(u(2));
dxdt = dxdt.';
end