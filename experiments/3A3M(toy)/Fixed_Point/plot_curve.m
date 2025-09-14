clc
close all
%%
figure('DefaultAxesFontSize',16)
avoid = ones(1,size(t,2));
avoid = avoid*59.1;
g = patch([0 t(end) t(end) 0],[58.8 58.8 59.1 59.1],'red', 'LineWidth', 2)
hold on
g.Annotation.LegendInformation.IconDisplayStyle = 'off';
target = ones(1,size(t,2));
target = target*(60-0.48);
h = patch([0 t(end) t(end) 0],[target(end) target(end) 60.4 60.4],'green', 'LineWidth', 2)
hold on
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(t,mac_spd(1,:)*60,'Color','blue', 'LineWidth', 2)
hold on
plot(t,mac_spd(2,:)*60,'Color','cyan', 'LineWidth', 2)
hold on
plot(t,mac_spd(3,:)*60,'Color','black', 'LineWidth', 2)
axis tight
%title('Frequency','Interpreter','latex')
ylabel('Frequency [Hz]','Interpreter','latex', 'FontSize',16)
xlabel('Time [s]','Interpreter','latex', 'FontSize',16)
%zlabel('Input','Interpreter','latex', 'FontSize',16)
leg1 = legend('Area 1','Area 2','Area 3');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',16);

figure('DefaultAxesFontSize',16)
plot(t,lmod_sig(1,:),'Color','blue', 'LineWidth', 2)
hold on
plot(t,lmod_sig(2,:),'Color','cyan', 'LineWidth', 2)
hold on
plot(t,lmod_sig(3,:),'Color','black', 'LineWidth', 2)
%axis tight
%title('Frequency','Interpreter','latex')
ylabel('$\Delta P$ [pu]','Interpreter','latex', 'FontSize',16)
xlabel('Time [s]','Interpreter','latex', 'FontSize',16)
%zlabel('Input','Interpreter','latex', 'FontSize',14)
leg1 = legend('Area 1','Area 2','Area 3');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',16);
