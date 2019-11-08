clearvars 
close all
fnum = 0;


tag = 'Test3Aa';
dat = importdata('./user102/Mag_v_t.txt');


sym = 4;  
tm = dat(:,1);
I0 = dat(:,2);
BHT0 = dat(:,3);
temHT0 = dat(:,4);
tauHT0 = dat(:,5);
rhocuHT0 = dat(:,6);
rhoHT0 = dat(:,7);

BHT0_cern = importdata('..\..\CERN_res\3Aa\B_avg_HT0_3Aa.csv');
rhocuHT0cern = importdata('..\..\CERN_res\3Aa\rhocu_avg_HT0_3Aa.csv');
rhoHT0cern = importdata('..\..\CERN_res\3Aa\rho_avg_HT0_3Aa.csv');
temHT0cern = importdata('..\..\CERN_res\3Aa\T_avg_HT0_3Aa.csv');
tauHT0cern = importdata('..\..\CERN_res\3Aa\tau_avg_HT0_3Aa.csv');



fnum=fnum+1;
h(fnum)=figure;
scatter(tm,BHT0,'filled')
hold on 
scatter(BHT0_cern(:,1),BHT0_cern(:,2),10,'filled');
% scatter(sQtot(:,1),sQtot(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Bave HT0 (T)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['Bave_HT0_vs_t_',tag],'-r300')
hold off



fnum=fnum+1;
h(fnum)=figure;
scatter(tm,temHT0,'filled')
hold on 
scatter(temHT0cern(:,1),temHT0cern (:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('temp ave HT0 (K)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['tempave_HT0_vs_t_',tag],'-r300')
hold off

fnum=fnum+1;
h(fnum)=figure;
scatter(tm,tauHT0,'filled')
hold on 
scatter(tauHT0cern(:,1),tauHT0cern(:,2)/2.0,10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('tau-IFCC ave HT0 (s)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['Tauave_HT0_vs_t_',tag],'-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(tm,rhocuHT0,'filled')
hold on 
scatter(rhocuHT0cern(:,1),rhocuHT0cern(:,2),10,'filled')
box on
hold on
grid on
title(tag)
% axis([0 1 1e-9 1e-7])
xlabel('tm (s)','FontSize',18)
ylabel('Cu rho ave HT0 (ohm-m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['rhocuave_HT0_vs_t_',tag],'-r300')
hold off



fnum=fnum+1;
h(fnum)=figure;
scatter(tm,rhoHT0,'filled')
hold on 
scatter(rhoHT0cern(:,1),rhoHT0cern(:,2),10,'filled')
box on
hold on
grid on
title(tag)
% axis([0 1 1e-9 1e-7])
xlabel('tm (s)','FontSize',18)
ylabel('rho ave HT0 (ohm-m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['rhoave_HT0_vs_t_',tag],'-r300')
hold off




% 
dlmwrite('ANSYS_B_avg_HT0_3Aa.csv',[tm,BHT0])
dlmwrite('ANSYS_temp_avg_HT0_3Aa.csv',[tm,temHT0])
dlmwrite('ANSYS_tau_avg_HT0_3Aa.csv',[tm,tauHT0])
dlmwrite('ANSYS_rhocu_avg_HT0_3Aa.csv',[tm,rhocuHT0])
dlmwrite('ANSYS_rho_avg_HT0_3Aa.csv',[tm,rhoHT0])

% dlmwrite('ANSYS_M_ifccy_avg_HT0_3Aa.csv',[tm Myave])

return 





