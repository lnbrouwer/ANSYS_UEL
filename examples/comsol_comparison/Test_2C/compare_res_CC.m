clearvars 
close all
fnum = 0;

%convert circu res
run .\user102_CC\write_trans_elem_res.m
datCirc = importdata('./user102_CC/circu_res.txt');
Icirc = datCirc(:,2);


tag = 'Test2CCC';

dat = importdata('./user102_CC/Mag_v_t.txt');
dat2 = importdata('./user102_CC/Therm_v_t.txt');
tmT = dat2(:,1);
TaveHT0 = dat2(:,2);


sym = 4;  
I0 = Icirc;
flux = dat(:,2)*sym;
LL = dat(:,3)*sym./I0;    %now divide by current here since didn't have before
BaveHT0 = dat(:,4);
tm = dat(:,5);
Mxave = dat(:,6);
Myave = dat(:,7);
JtauHT0 = dat(:,8);
JtauALL = dat(:,9)*sym;
JsumHT0 = dat(:,10);
JALL = dat(:,11)*sym;
rhoHT0 = dat(:,12);








sQHT0 = importdata('..\CERN_res\2C\Q_avg_HT0_2C2.csv');
sQtot = importdata('..\CERN_res\2C\Q_tot_2C2.csv');
sRhoHT0 = importdata('..\CERN_res\2C\rho_avg_HT0_2C2.csv');
sTHT0 = importdata('..\CERN_res\2C\T_avg_HT0_2C2.csv');


fnum=fnum+1;
h(fnum)=figure;
scatter(tm,JsumHT0,'filled')
hold on 
scatter(sQHT0(:,1),sQHT0(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Qtot HT0 (W/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['Qtot_HT0_vs_t_',tag],'-r300')
hold off

fnum=fnum+1;
h(fnum)=figure;
scatter(tm,JALL,'filled')
hold on 
scatter(sQtot(:,1),sQtot(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Qtot all (W/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['Qtot_all_vs_t_',tag],'-r300')
hold off

fnum=fnum+1;
h(fnum)=figure;
scatter(tm,rhoHT0,'filled')
hold on 
scatter(sRhoHT0(:,1),sRhoHT0(:,2),10,'filled')
% axis([sRhoHT0(1,1) sRhoHT0(length(sRhoHT0(:,1)),1) 0.95*min(sRhoHT0(:,2)) 1.05*max(sRhoHT0(:,2))])
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('ave rho HT0 (ohm-m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['rho_HT0_vs_t_',tag],'-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(tmT,TaveHT0,'filled')
hold on 
scatter(sTHT0(:,1),sTHT0(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('ave temp HT0 (K)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['temp_HT0_vs_t_',tag],'-r300')
hold off





dlmwrite('ANSYS_Q_avg_HT0_2C_CC.csv',[tm JsumHT0])
dlmwrite('ANSYS_Q_tot_2C_CC.csv',[tm JALL])
dlmwrite('ANSYS_rho_avg_HT0_2C_CC.csv',[tm rhoHT0])
dlmwrite('ANSYS_T_avg_HT0_2C_CC.csv',[tm TaveHT0])



return 






