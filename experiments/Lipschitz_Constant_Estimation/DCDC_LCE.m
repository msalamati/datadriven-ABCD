%author: B.Wooding
%date: Jan 2022
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
dim = 2; % state dimensions
u = [0; 1]; % list of all inputs
ins = size(u,1); % number of inputs
M = zeros(ins,1); % store lipschitz constant estimations for each input
lb = [0.65 4.95]; %lowerbound of system domain
ub = [1.65 5.95]; %upperbound of system domain

% uniform distributions [a,b] x [a,b] for x_i and y_i
dist = {
    makedist('Uniform','lower',lb(1),'upper',ub(1));
    makedist('Uniform','lower',lb(2),'upper',ub(2));
    };
%% Implementation
% loop through all inputs
for ix = 1:ins
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
            %calculate s_i
            [t1 ode1]=ode45(@g,[0 .5], xx, [],u(ix));
            [t2 ode2]=ode45(@g,[0 .5], yy, [],u(ix));
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
    M(ix) = params(3); % save location parameter to set of Lipschitz Constants
end
%% Output
max(M)
%% functions
function dxdt = g(t,x,u)
h  = 10*2/4e3;
r0 = 1.0;
vs = 1.0;
rl = 0.05;
rc = rl / 10;
xl = 3.0;
xc = 70.0;
B = [vs/xl;0];

A = zeros(2);
if(u==1)
    A(1,1) = -rl / xl;
    A(1,2) = 0;
    A(2,1) = 0;
    A(2,2) = (-1 / xc) * (1 / (r0 + rc));
else
    A(1,1) = (-1 / xl) * (rl + ((r0 * rc) / (r0 + rc)));
    A(1,2) =  ((-1 / xl) * (r0 / (r0 + rc))) / 5;
    A(2,1) = 5 * (r0 / (r0 + rc)) * (1 / xc);
    A(2,2) = (-1 / xc) * (1 / (r0 + rc));
end
dxdt = A*x + B;
end