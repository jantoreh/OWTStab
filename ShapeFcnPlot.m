figure(1);
subplot(1,2,1);
x=[0:0.5:50];
plot(x,BendingModeFun(x,1),'-k');
set(gca,'FontSize',10);
xlabel('Radius, $z$ [m]','Interpreter','latex','fontsize',12);
ylabel('Bending shape fcn, $\phi_b$ [m]','Interpreter','latex','fontsize',12);
xlim([0 50]);
%set(gca,'YTick',-0.5:0.5:0.5)
%set(gca,'XTick',-60:10:120)
%set(gca,'YGrid','on');
set(gca,'PlotBoxAspectRatio',[2 1 1])
subplot(1,2,2);
plot(x,TorsionalModeFun(x,1),'-k');
set(gca,'FontSize',10);
xlabel('Radius, $z$ [m]','Interpreter','latex','fontsize',12);
ylabel('Torsional shape fcn, $\phi_t$ [rad]','Interpreter','latex','fontsize',12);
%ylim([0 5]);
xlim([0 50]);
%set(gca,'YTick',0:1:5)
%set(gca,'XTick',-60:10:120)
set(gca,'PlotBoxAspectRatio',[2 1 1])
set(gcf,'Units','centimeters')
set(gcf,'Position',[1,1,26,12])
set(gcf,'PaperPositionMode','auto')
saveas(gcf,'shapefcns','eps')
