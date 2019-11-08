clearvars 
close all
fnum = 0;

%convert circu res
run .\user102_CC\write_trans_elem_res.m
datCirc = importdata('./user102_CC/circu_res.txt');
Icirc = datCirc(:,2);


tag = 'Test3AcCC';
dat = importdata('./user102_CC/Mag_v_t.txt');

sym = 4;  
tm = dat(:,1);
I0 = Icirc;
BHT0 = dat(:,3);
temHT0 = dat(:,4);
tauHT0 = dat(:,5);
rhocuHT0 = dat(:,6);
rhoHT0 = dat(:,7);


BHT0_cern = importdata('..\..\CERN_res\3Ac\B_avg_HT0_3Ac.csv');
rhocuHT0_cern = importdata('..\..\CERN_res\3Ac\rhocu_avg_HT0_3Ac.csv');
rhoHT0_cern = importdata('..\..\CERN_res\3Ac\rho_avg_HT0_3Ac.csv');
temHT0_cern = importdata('..\..\CERN_res\3Ac\T_avg_HT0_3Ac.csv');
tauHT0_cern = importdata('..\..\CERN_res\3Ac\tau_avg_HT0_3Ac.csv');



fnum=fnum+1;
h(fnum)=figure;
scatter(tm,BHT0,'filled')
hold on 
scatter(BHT0_cern(:,1),BHT0_cern(:,2),10,'filled')
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
scatter(temHT0_cern(:,1),temHT0_cern(:,2),10,'filled')
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
scatter(tauHT0_cern(:,1),tauHT0_cern(:,2)/2.0,10,'filled')
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
scatter(rhocuHT0_cern(:,1),rhocuHT0_cern(:,2),10,'filled')
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
scatter(rhoHT0_cern(:,1),rhoHT0_cern(:,2),10,'filled')
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
dlmwrite('ANSYS_B_avg_HT0_3Ac_CC.csv',[tm,BHT0])
dlmwrite('ANSYS_temp_avg_HT0_3Ac_CC.csv',[tm,temHT0])
dlmwrite('ANSYS_tau_avg_HT0_3Ac_CC.csv',[tm,tauHT0])
dlmwrite('ANSYS_rhocu_avg_HT0_3Ac_CC.csv',[tm,rhocuHT0])
dlmwrite('ANSYS_rho_avg_HT0_3Ac_CC.csv',[tm,rhoHT0])

% dlmwrite('ANSYS_M_ifccy_avg_HT0_3Aa.csv',[tm Myave])

return 





