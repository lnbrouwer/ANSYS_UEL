clearvars 
close all
fnum = 0;


tag = 'SW01';
dat = importdata('./user102/Mag_v_t.txt');


mu = 4*pi*1e-7;

dBdt = dat(1,1);
tau = dat(1,2);
fcond = dat(1,3);
HaveHT0 = dat(:,4);
BaveHT0 = dat(:,5);
tm = dat(:,6);
Mxave = dat(:,7);
Myave = dat(:,8);
Pave = dat(:,9);



Bi = dBdt*(tm + tau*(exp(-tm/tau)-1));
Mana = -2*dBdt*tau*(1-exp(-tm/tau))/mu;
Hana = dBdt*(tm + tau*(1-exp(-tm/tau)))/mu;
Mmax = -2*tau*dBdt/mu;
Pana = 2*dBdt*dBdt*tau*(1-exp(-tm/tau)).*(1-exp(-tm/tau))/mu;

fnum=fnum+1;
h(fnum)=figure;
scatter(tm*1000,BaveHT0,'filled')
hold on 
plot(tm*1000,Bi,'linewidth',1.5)
box on
hold on
% grid on
title(tag)
xlabel('time (ms)','FontSize',18)
ylabel('Bi (T)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
% legend('ANSYS','Analytic','Location','NorthEastOutside')
legend('ANSYS','Eqn. 9','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['Bi_vs_t_',tag],'-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(tm*1000,Myave,'filled')
hold on
plot(tm*1000,Mana,'linewidth',1.5)
% plot([tm(1),tm(end)],[Mmax,Mmax],'--k')
box on
hold on
% grid on
title(tag)
xlabel('time (ms)','FontSize',18)
ylabel('Me (A/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
% legend('ANSYS','Analytic','Location','NorthEastOutside')
legend('ANSYS','Eqn. 10','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['Me_vs_t_',tag],'-r300')
hold off

fnum=fnum+1;
h(fnum)=figure;
scatter(tm*1000,HaveHT0,'filled')
hold on 
plot(tm*1000,Hana,'linewidth',1.5)
box on
hold on
% grid on
title(tag)
xlabel('time (ms)','FontSize',18)
ylabel('Hi (A/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
% legend('ANSYS','Analytic','Location','NorthEastOutside')
legend('ANSYS','Eqn. 11','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['Hi_vs_t_',tag],'-r300')
hold off




fnum=fnum+1;
h(fnum)=figure;
scatter(tm*1000,Pave,'filled')
hold on
plot(tm*1000,Pana,'linewidth',1.5)
% plot([tm(1),tm(end)],[Mmax,Mmax],'--k')
box on
hold on
% grid on
title(tag)
xlabel('time (ms)','FontSize',18)
ylabel('Pe (W/m^3)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
% legend('ANSYS','Analytic','Location','NorthEastOutside')
legend('ANSYS','Eqn. 12','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['Pe_vs_t_',tag],'-r300')
hold off










