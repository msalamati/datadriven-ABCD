%
% vehicle.m
%
% created on: 09.10.2015
%     author: rungger
% editted on 14.09.2021
%     edited by: wooding
%
% see readme file for more information on the vehicle example
%
% you need to run ./vehicleModel_c++ and ./vehicleDDS_c++ binaries first 
%
% so that the files: vehicle_ss_DDS.bdd 
%                    vehicle_obst_DDS.bdd
%                    vehicle_target_DDS.bdd
%                    vehicle_controller_DDS.bdd
%                    vehicle_controller_Model.bdd
% are created
%

function vehicle
clear set
close all



%% simulation
% target set
lb=[9 0];
ub=lb+0.5;
v=[9 0; 9.51  0; 9 0.51; 9.51 0.51];
% initial state
x0=[0, 1.2, 0];

%data driven controller printed in black
controller=SymbolicSet('Controllers/vehicle_controller_DDS.bdd','projection',[1 2 3]);
%model based controller printed in red
controller2=SymbolicSet('Controllers/vehicle_controller_Model.bdd','projection',[1 2 3]);
%target region
target=SymbolicSet('Target_Obst_SS/vehicle_target_DDS.bdd');

y=x0;
y2=x0;
v=[];
v2=[];

%plot path for data driven control vehicle
while(1)
  
  if (target.isElement(y(end,:)))
    break;
  end 

  u=controller.getInputs(y(end,:));
  v=[v; u(1,:)];
  [t x]=ode45(@vehicle_ode,[0 .3], y(end,:), [],u(1,:));
  x(end,:);
  y=[y; x(end,:)];
end
disp('DD complete')
%plot path for model based control vehicle
while(1)
  
  if (target.isElement(y2(end,:)))
    break;
  end 
  
  u2= controller2.getInputs(y2(end,:));
  v2=[v2; u2(1,:)];
  [t2 x2]=ode45(@vehicle_ode,[0 .3], y2(end,:), [],u2(1,:));
  x2(end,:);
  y2=[y2; x2(end,:)];
end
disp('MB complete')
%% plot the vehicle domain
% colors
colors=get(groot,'DefaultAxesColorOrder');
figure('DefaultAxesFontSize',16)

% load the symbolic set containig the abstract state space
set=SymbolicSet('Target_Obst_SS/vehicle_ss_DDS.bdd','projection',[1 2]);
plotCells(set,'facecolor','none','edgec',[0.8 0.8 0.8],'linew',.1)
hold on

% load the symbolic set containig obstacles
set=SymbolicSet('Target_Obst_SS/vehicle_obst_DDS.bdd','projection',[1 2]);
plotCells(set,'facecolor',colors(1,:)*0.5+0.5,'edgec',colors(1,:),'linew',.1)

% load the symbolic set containing target set
set=SymbolicSet('Target_Obst_SS/vehicle_target_DDS.bdd','projection',[1 2]);
plotCells(set,'facecolor',colors(2,:)*0.5+0.5,'edgec',colors(2,:),'linew',.1)

% plot initial state  and trajectory
plot(y(:,1),y(:,2),'k.-','LineWidth',1.5) %plot data driven controller
plot(y2(:,1),y2(:,2),'r.-','LineWidth',1.5) %plot model based controller
plot(y(1,1),y(1,2),'.','color',colors(5,:),'markersize',20)


box on
axis([-.5 10.5 -.5 10.5])
axis tight
ylabel('y(m)','Interpreter','latex', 'FontSize',18)
xlabel('x(m)','Interpreter','latex', 'FontSize',18)

end

function dxdt = vehicle_ode(t,x,u)

  dxdt = zeros(3,1);
  c=atan(tan(u(2))/2);

  dxdt(1)=u(1)*cos(c+x(3))/cos(c);
  dxdt(2)=u(1)*sin(c+x(3))/cos(c);
  dxdt(3)=u(1)*tan(u(2));


end

