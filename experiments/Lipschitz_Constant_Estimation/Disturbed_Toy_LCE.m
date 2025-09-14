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

rng(1); % set random seed to 1
%% Parameters
n = 3; % number of slopes
m = 5; % number of data points for RWD
dim = 3; % state dimensions
u = [0; 0.025; 0.05; 0.075; 0.1; 0.125; 0.15; 0.175; 0.2; 0.225; 0.25; 0.275; 0.3; 0.325; 0.35; 0.375; 0.4; 0.425; 0.45; 0.475; 0.5]; % list of all inputs
ins = size(u,1); % number of inputs
M = zeros(ins,1); % store lipschitz constant estimations for each input
eta = 0.005; %0.0015;
dim1 = [-0.02:eta:0.02];
dim2 = [-0.05:eta:0.05];
[p,q] = meshgrid(dim1, dim2);
pairs = [p(:) q(:)];
LSol = zeros(size(pairs,1),1);
%lb = [-0.02 -0.05 -0.12]; %lowerbound of system domain
%ub = [ 0.02  0.05  0.12]; %upperbound of system domain
% uniform distributions [a,b] x [a,b] for x_i and y_i

for pair = 1:size(pairs,1)
dist = {
        makedist('Uniform','lower',pairs(pair,1),'upper',pairs(pair,1)+eta);
        makedist('Uniform','lower',pairs(pair,2),'upper',pairs(pair,2)+eta);
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
            for i=1:dim-1
                % Choose random x_i and y_i
                xx(i,:)=random(dist{i});
                yy(i,:)=random(dist{i});
            end
            
            [t1 ode1]=ode45(@g,[0 .5], xx, [],u(ix));
            [t2 ode2]=ode45(@g,[0 .5], yy, [],u(ix));
                %g(xx,u(ix))
                %g(yy,u(ix))
                slopes(j) = norm(ode1(end,:) - ode2(end,:))/norm(xx-yy);

            %calculate s_i
            %slopes(j) = norm(g(xx,u(ix)) - g(yy,u(ix)))/norm(xx-yy);
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
    %0.76482 5.61034 max(l)
    M(ix) = params(3); % save location parameter to set of Lipschitz Constants
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
%zlim([0.9 1.1])

%% functions
function dxdt = g(t,x,u)
A = [0.00027563  0       0;
    0          -0.3951  0.687;
    0          -0.6869 -0.016];

B = [0.00031166; 0.1359; 0.0230];

E = [0.00033103 0.00031244;
     0.1309     0.1308;
     0.0250     0.0233];
 
 w = [0.3; 0.2];

dxdt = A*x + B*u + E*w;
end