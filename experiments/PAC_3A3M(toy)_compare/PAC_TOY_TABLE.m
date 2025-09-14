clc
clear all
%close all
%% Compare Models
D_d = SymbolicSet('Disturbed_Controllers/Toy.bdd'); %dcdc disturbed
%D_u = SymbolicSet('Undisturbed_Controllers/PAC.bdd');%dcdc undisturbed

P_d = SymbolicSet('Disturbed_Controllers/PAC.bdd');%PAC disturbed
%P_u = SymbolicSet('Undisturbed_Controllers/PAC_controller.bdd');%PAC undisturbed

M_d = SymbolicSet('Disturbed_Controllers/Model.bdd');%model_disturbed
%M_u = SymbolicSet('Undisturbed_Controllers/dcdcModel_controller.bdd');%model undisturbed

% PAC points
PAC_d = P_d.points();
PAC_d = unique(PAC_d(:,1:3),'rows');
%size(PAC_d)
%PAC_u = P_u.points();
% RWD points
RWD_d = D_d.points();
RWD_d = unique(RWD_d(:,1:3),'rows');
%size(RWD_d)
%RWD_u = D_u.points();
% Model points
Model_d = M_d.points();
Model_d = unique(Model_d(:,1:3),'rows');
%size(Model_d)
%Model_u = M_u.points();
%% Calc Set Diff
RWD_Model = setdiff(RWD_d,Model_d,'rows');
%size(RWD_Model)/size(RWD_d)*100
RWD_PAC = setdiff(RWD_d,PAC_d,'rows');
%size(RWD_PAC)/size(RWD_d)*100
PAC_Model = setdiff(PAC_d,Model_d,'rows');
%size(PAC_Model)/size(PAC_d)*100
PAC_RWD = setdiff(PAC_d,RWD_d,'rows');
%size(PAC_RWD)/size(PAC_d)*100
Model_PAC = setdiff(Model_d,PAC_d,'rows');
%size(Model_PAC)/size(Model_d)*100
Model_RWD = setdiff(Model_d,RWD_d,'rows');
%size(Model_RWD)/size(Model_d)*100
%% Calc Intercept
I_Model_RWD = intersect(RWD_d,Model_d,'rows');
%size(I_Model_RWD)/size(Model_d)*100
I_RWD_PAC = intersect(RWD_d,PAC_d,'rows');
%size(I_RWD_PAC)/size(PAC_d)*100
I_Model_PAC = intersect(PAC_d,Model_d,'rows');
%size(I_Model_PAC)/size(Model_d)*100

%% For Loop exception PAC
a = 0;
b = 0;
for i = 1:size(PAC_Model,1)
    xb = PAC_Model(i,:) +[0.0014/2 0.0014/2 0.0014/2];
    z = P_d.getInputs(xb);
    for j = 1:size(P_d.getInputs(xb))
    u = z(j);
    [t1 x1]= ode45(@g, [0 0.4], xb, [], u);
    out=x1(end,:);
    %v = min(P_d.getInputs(out));
    %[ i out v]
    %[t1 x1]= ode45(@g, [0 0.4], out, [], v);
    %out=x1(end,:);
    %[u xb out]
    %w = min(P_d.getInputs(out));
    try
        w = min(P_d.getInputs(out));
        %[out w]
    catch
        %warning('Problem using function.  Assigning a value of 0.');
        a = a + 1;
    end
    b = b + 1;
    end
end
%% function
function dxdt = g(t,x,u)
A = [0.00027563  0       0;
    0          -0.3951  0.687;
    0          -0.6869 -0.016];

B = [0.00031166; 0.1359; 0.0230];

E = [0.00033103 0.00031244;
    0.1309     0.1308;
    0.0250     0.0233];

w = [0.3; 0.2];

dxdt = A*x+ B*-u + E*w;
end