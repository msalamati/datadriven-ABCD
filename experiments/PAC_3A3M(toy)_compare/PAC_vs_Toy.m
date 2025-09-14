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
PAC_d = PAC_d(:,1:3);
%PAC_u = P_u.points();
% RWD points
RWD_d = D_d.points();
RWD_d = RWD_d(:,1:3);
%RWD_u = D_u.points();
% Model points
Model_d = M_d.points();
Model_d = Model_d(:,1:3);
%Model_u = M_u.points();
%% Calc Set Diff
RWD_Model = setdiff(RWD_d,Model_d,'rows');
%size(RWD_Model)
RWD_PAC = setdiff(RWD_d,PAC_d,'rows');
%size(RWD_PAC)
PAC_Model = setdiff(PAC_d,Model_d,'rows');
%size(PAC_Model)
PAC_RWD = setdiff(PAC_d,RWD_d,'rows');
%size(PAC_RWD)
Model_PAC = setdiff(Model_d,PAC_d,'rows');
%size(Model_PAC)
Model_RWD = setdiff(Model_d,RWD_d,'rows');
%size(Model_RWD)
%% Calc Intercept
I_Model_RWD = intersect(RWD_d,Model_d,'rows');
%size(I_RWD_Model)
I_RWD_PAC = intersect(RWD_d,PAC_d,'rows');
%size(I_RWD_PAC)
I_Model_PAC = intersect(PAC_d,Model_d,'rows');
%size(I__Model_PAC)

%% Calc Union
% U_RWD_Model = union(RWD_d,Model_d,'rows');
% %size(U_RWD_Model)
% U_RWD_PAC = union(RWD_d,PAC_d,'rows');
% %size(U_RWD_PAC)
% U_PAC_Model = union(PAC_d,Model_d,'rows');
%size(U_PAC_Model)

%% Comparison of Set Difference 2d
figure('DefaultAxesFontSize',16)
plot(Model_d(:,1),Model_d(:,2),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on
plot(PAC_ds(:,1),PAC_ds(:,2),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','red','MarkerFaceColor','red')
hold on
plot(RWD_ds(:,1),RWD_ds(:,2),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','black','MarkerFaceColor','black')
axis tight
%title('Disturbed DCDC','Interpreter','latex')
ylabel('Dim 2','Interpreter','latex', 'FontSize',18)
xlabel('Dim 1','Interpreter','latex', 'FontSize',18)
%zlabel('Dim 3','Interpreter','latex', 'FontSize',18)
leg1 = legend('Model','PAC','RSA');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',16);
set(gca,'xtick',[])
set(gca,'ytick',[])

figure('DefaultAxesFontSize',16)
plot(Model_d(:,1),Model_d(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on
plot(PAC_ds(:,1),PAC_ds(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','red','MarkerFaceColor','red')
hold on
plot(RWD_ds(:,1),RWD_ds(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','black','MarkerFaceColor','black')
axis tight
%title('Disturbed DCDC','Interpreter','latex')
ylabel('Dim 3','Interpreter','latex', 'FontSize',18)
xlabel('Dim 1','Interpreter','latex', 'FontSize',18)
%zlabel('Dim 3','Interpreter','latex', 'FontSize',18)
leg1 = legend('Model','PAC','RSA');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',16);
set(gca,'xtick',[])
set(gca,'ytick',[])

figure('DefaultAxesFontSize',16)
plot(Model_d(:,2),Model_d(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on
plot(PAC_ds(:,2),PAC_ds(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','red','MarkerFaceColor','red')
hold on
plot(RWD_ds(:,2),RWD_ds(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','black','MarkerFaceColor','black')
axis tight
%title('Disturbed DCDC','Interpreter','latex')
ylabel('Dim 3','Interpreter','latex', 'FontSize',18)
xlabel('Dim 2','Interpreter','latex', 'FontSize',18)
%zlabel('Dim 3','Interpreter','latex', 'FontSize',18)
leg1 = legend('Model','PAC','RSA');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',16);
set(gca,'xtick',[])
set(gca,'ytick',[])

%% Comparison of PAC and RWD
figure('DefaultAxesFontSize',16)
plot(I_RWD_PAC(:,1),I_RWD_PAC(:,2),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on
plot(PAC_RWD(:,1),PAC_RWD(:,2),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','red','MarkerFaceColor','red')
hold on
plot(RWD_PAC(:,1),RWD_PAC(:,2),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','black','MarkerFaceColor','black')
axis tight
%title('Disturbed DCDC','Interpreter','latex')
ylabel('Dim 2','Interpreter','latex', 'FontSize',18)
xlabel('Dim 1','Interpreter','latex', 'FontSize',18)
%zlabel('Dim 3','Interpreter','latex', 'FontSize',18)
leg1 = legend('PAC $\cap$ RSA','PAC $\backslash$ RSA','RSA $\backslash$ PAC');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',16);
set(gca,'xtick',[])
set(gca,'ytick',[])

figure('DefaultAxesFontSize',16)
plot(I_RWD_PAC(:,1),I_RWD_PAC(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on
plot(PAC_RWD(:,1),PAC_RWD(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','red','MarkerFaceColor','red')
hold on
plot(RWD_PAC(:,1),RWD_PAC(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','black','MarkerFaceColor','black')
axis tight
%title('Disturbed DCDC','Interpreter','latex')
ylabel('Dim 3','Interpreter','latex', 'FontSize',18)
xlabel('Dim 1','Interpreter','latex', 'FontSize',18)
%zlabel('Dim 3','Interpreter','latex', 'FontSize',18)
leg1 = legend('PAC $\cap$ RSA','PAC $\backslash$ RSA','RSA $\backslash$ PAC');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',16);
set(gca,'xtick',[])
set(gca,'ytick',[])

figure('DefaultAxesFontSize',16)
plot(I_RWD_PAC(:,2),I_RWD_PAC(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on
plot(PAC_RWD(:,2),PAC_RWD(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','red','MarkerFaceColor','red')
hold on
plot(RWD_PAC(:,2),RWD_PAC(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','black','MarkerFaceColor','black')
axis tight
%title('Disturbed DCDC','Interpreter','latex')
ylabel('Dim 3','Interpreter','latex', 'FontSize',18)
xlabel('Dim 2','Interpreter','latex', 'FontSize',18)
%zlabel('Dim 3','Interpreter','latex', 'FontSize',18)
leg1 = legend('PAC $\cap$ RSA','PAC $\backslash$ RSA','RSA $\backslash$ PAC');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',16);
set(gca,'xtick',[])
set(gca,'ytick',[])

%% Comparison of PAC and Model
figure('DefaultAxesFontSize',16)
plot(I_Model_PAC(:,1),I_Model_PAC(:,2),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on
plot(PAC_Model(:,1),PAC_Model(:,2),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','red','MarkerFaceColor','red')
hold on
plot(Model_PAC(:,1),Model_PAC(:,2),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','black','MarkerFaceColor','black')
%hold on
%plot(0.0003,0.0611,'.','color','green','markersize',20)
axis tight
%title('Disturbed DCDC','Interpreter','latex')
ylabel('Dim 2','Interpreter','latex', 'FontSize',18)
xlabel('Dim 1','Interpreter','latex', 'FontSize',18)
%zlabel('Dim 3','Interpreter','latex', 'FontSize',18)
leg1 = legend('PAC $\cap$ Model','PAC $\backslash$ Model','Model $\backslash$ PAC');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',16);
set(gca,'xtick',[])
set(gca,'ytick',[])

figure('DefaultAxesFontSize',16)
plot(I_Model_PAC(:,1),I_Model_PAC(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on
plot(PAC_Model(:,1),PAC_Model(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','red','MarkerFaceColor','red')
hold on
plot(Model_PAC(:,1),Model_PAC(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','black','MarkerFaceColor','black')
%hold on
%plot(0.0003,-0.0150,'.','color','green','markersize',20)
axis tight
%title('Disturbed DCDC','Interpreter','latex')
ylabel('Dim 3','Interpreter','latex', 'FontSize',18)
xlabel('Dim 1','Interpreter','latex', 'FontSize',18)
%zlabel('Dim 3','Interpreter','latex', 'FontSize',18)
leg1 = legend('PAC $\cap$ Model','PAC $\backslash$ Model','Model $\backslash$ PAC');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',16);
set(gca,'xtick',[])
set(gca,'ytick',[])

figure('DefaultAxesFontSize',16)
plot(I_Model_PAC(:,2),I_Model_PAC(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on
plot(PAC_Model(:,2),PAC_Model(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','red','MarkerFaceColor','red')
hold on
plot(Model_PAC(:,2),Model_PAC(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','black','MarkerFaceColor','black')
%hold on
%plot(0.0611,-0.015,'.','color','green','markersize',20)
axis tight
%title('Disturbed DCDC','Interpreter','latex')
ylabel('Dim 3','Interpreter','latex', 'FontSize',18)
xlabel('Dim 2','Interpreter','latex', 'FontSize',18)
%zlabel('Dim 3','Interpreter','latex', 'FontSize',18)
leg1 = legend('PAC $\cap$ Model','PAC $\backslash$ Model','Model $\backslash$ PAC');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',16);
set(gca,'xtick',[])
set(gca,'ytick',[])

%% Comparison of RWD and Model
figure('DefaultAxesFontSize',16)
plot(I_Model_RWD(:,1),I_Model_RWD(:,2),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on
plot(RWD_Model(:,1),RWD_Model(:,2),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','red','MarkerFaceColor','red')
hold on
plot(Model_RWD(:,1),Model_RWD(:,2),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','black','MarkerFaceColor','black')
axis tight
%title('Disturbed DCDC','Interpreter','latex')
ylabel('Dim 2','Interpreter','latex', 'FontSize',18)
xlabel('Dim 1','Interpreter','latex', 'FontSize',18)
%zlabel('Dim 3','Interpreter','latex', 'FontSize',18)
leg1 = legend('RSA $\cap$ Model','RSA $\backslash$ Model','Model $\backslash$ RSA');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',16);
set(gca,'xtick',[])
set(gca,'ytick',[])

figure('DefaultAxesFontSize',16)
plot(I_Model_RWD(:,1),I_Model_RWD(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on
plot(RWD_Model(:,1),RWD_Model(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','red','MarkerFaceColor','red')
hold on
plot(Model_RWD(:,1),Model_RWD(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','black','MarkerFaceColor','black')
axis tight
%title('Disturbed DCDC','Interpreter','latex')
ylabel('Dim 3','Interpreter','latex', 'FontSize',18)
xlabel('Dim 1','Interpreter','latex', 'FontSize',18)
%zlabel('Dim 3','Interpreter','latex', 'FontSize',18)
leg1 = legend('RSA $\cap$ Model','RSA $\backslash$ Model','Model $\backslash$ RSA');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',16);
set(gca,'xtick',[])
set(gca,'ytick',[])

figure('DefaultAxesFontSize',16)
plot(I_Model_RWD(:,2),I_Model_RWD(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on
plot(RWD_Model(:,2),RWD_Model(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','red','MarkerFaceColor','red')
hold on
plot(Model_RWD(:,2),Model_RWD(:,3),'.','Marker','square','MarkerSize',15,'MarkerEdgeColor','black','MarkerFaceColor','black')
axis tight
%title('Disturbed DCDC','Interpreter','latex')
ylabel('Dim 3','Interpreter','latex', 'FontSize',18)
xlabel('Dim 2','Interpreter','latex', 'FontSize',18)
%zlabel('Dim 3','Interpreter','latex', 'FontSize',18)
leg1 = legend('RSA $\cap$ Model','RSA $\backslash$ Model','Model $\backslash$ RSA');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',16);
set(gca,'xtick',[])
set(gca,'ytick',[])

%% 3D Model RSA
figure('DefaultAxesFontSize',16)
plot3(I_Model_RWD(:,1),I_Model_RWD(:,2),I_Model_RWD(:,3),'.')
hold on
plot3(RWD_Model(:,1),RWD_Model(:,2),RWD_Model(:,3),'.')
hold on
plot3(Model_RWD(:,1),Model_RWD(:,2),Model_RWD(:,3),'.')
axis tight
%title('Disturbed DCDC','Interpreter','latex')
ylabel('Dim 2','Interpreter','latex', 'FontSize',18)
xlabel('Dim 1','Interpreter','latex', 'FontSize',18)
zlabel('Dim 3','Interpreter','latex', 'FontSize',18)
leg1 = legend('RSA $\cap$ Model','RSA $\backslash$ Model','Model $\backslash$ RSA');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',16);

%% 3D Model PAC
figure('DefaultAxesFontSize',16)
plot3(I_Model_PAC(:,1),I_Model_PAC(:,2),I_Model_PAC(:,3),'.')
hold on
plot3(PAC_Model(:,1),PAC_Model(:,2),PAC_Model(:,3),'.')
hold on
plot3(Model_PAC(:,1),Model_PAC(:,2),Model_PAC(:,3),'.')
axis tight
%title('Disturbed DCDC','Interpreter','latex')
ylabel('Dim 2','Interpreter','latex', 'FontSize',18)
xlabel('Dim 1','Interpreter','latex', 'FontSize',18)
zlabel('Dim 3','Interpreter','latex', 'FontSize',18)
leg1 = legend('PAC $\cap$ Model','PAC $\backslash$ Model','Model $\backslash$ PAC');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',16);

%% 3D RSA PAC
figure('DefaultAxesFontSize',16)
plot3(I_RWD_PAC(:,1),I_RWD_PAC(:,2),I_RWD_PAC(:,3),'.')
hold on
plot3(PAC_RWD(:,1),PAC_RWD(:,2),PAC_RWD(:,3),'.')
hold on
plot3(RWD_PAC(:,1),RWD_PAC(:,2),RWD_PAC(:,3),'.')
axis tight
%title('Disturbed DCDC','Interpreter','latex')
ylabel('Dim 2','Interpreter','latex', 'FontSize',18)
xlabel('Dim 1','Interpreter','latex', 'FontSize',18)
zlabel('Dim 3','Interpreter','latex', 'FontSize',18)
leg1 = legend('PAC $\cap$ RSA','PAC $\backslash$ RSA','RSA $\backslash$ PAC');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',16);

