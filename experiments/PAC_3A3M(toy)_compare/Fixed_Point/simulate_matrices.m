%% Desired System Matrices
load sysFull;  % previously saved system
sw_con(1,7)=0.01;
q = 3;   % number of inputs
p = 1;   % number of outputs
n = 63; % state dimension
% sysFull = drss(n,p, q);
r = 3;  % reduced model order
%% Reduce Model Order
%sysb = sysFull;
 [sysb,g,Ti,T] = balreal(sysFull); % balance system
 % check: https://people.kth.se/~hsan/modred_files/Lecture5.pdf
 
 A = sysb.A([1:r],[1:r]);
 B = sysb.B([1:r],1);
 C = sysb.C(:,[1:r]);
 E = sysb.B([1:r],2:3);

%sysBT = modred(sysb,r+1:n,'truncate'); % Truncate balanced system
%% simulation values
h = sw_con(1,7); %timestep
T = 50; %simulation length
x = zeros(size(A,1),1);
u = [0]; %tg_sig (input)
e = [3;2]; %lmod_sig (disturbance)

[t, x1] = ode45(@sim, [0 T],x,[],A,B,u,E,e);
y = x1.'; % state values
figure
plot(t,y) %plot frequency over time (Hz)
figure
plot(t,C*y)
%% function
function dydt=sim(t,x,A,B,u,E,e)
if(t>1)
    dydt = A*x+B*u+E*e;
else
    dydt = A*x+B*u;
end
end

