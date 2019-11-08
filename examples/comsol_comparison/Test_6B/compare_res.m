clearvars 
close all
fnum = 0;



%convert circu res
run .\user102_CC\write_trans_elem_res.m
datCirc = importdata('./user102_CC/circu_res.txt');    
% dlmwrite('circu_res.txt',[tm I1 I2 IR IRC IS ISC V1 V2 VR VRC VS VSC]);
tm = datCirc(:,1);

I1 = datCirc(:,2);
I2 = datCirc(:,3);
IR = datCirc(:,4);
IRC = datCirc(:,5);
IS = datCirc(:,6);
ISC = datCirc(:,7);


V1 = datCirc(:,8);
V2 = datCirc(:,9);
VR = datCirc(:,10);
VRC = datCirc(:,11);
VS = datCirc(:,12);
VSC = datCirc(:,13);
VC = datCirc(:,14);
IC = datCirc(:,15);



dat = importdata('./user102_CC/Mag_v_t.txt');
dat2 = importdata('./user102_CC/Therm_v_t.txt');
tmT = dat2(:,1);
TaveHT0 = dat2(:,2);


dat3 = importdata('./user102_CC/Thot_v_t.txt');
tmTH = dat3(:,1);
Thot = dat3(:,2);

tag = 'Test6B';

% 
sym = 4;  
Li = 9.20
Lc = 10.11


EE = dat(:,3);
Ei = (EE(1)+dat(1,4))*sym*Li;
Ef = (EE(length(EE))+dat(1,5))*sym*Li;
dE = Ei-Ef;
Emiss = Ef/Ei;


tauHT0 = dat(:,6);
MxHT0 = dat(:,7);
MyHT0 = dat(:,8);
PtauHT0 = dat(:,9);
PtauSUM = dat(:,10);

PresHT0 = dat(:,11);
PresSUM = dat(:,12);
Rcoil = dat(:,13)*sym;
IFCU = dat(:,14);


% 
% %integrate the dump power to find the loss
% pp = interp1(tm,abs(PR),'linear','pp');
% f1 = @(x) ppval(pp,x);
% QRtot = quad(f1,5e-3,tm(length(tm)))
% 
% 
% xint = 5e-3:5e-3:tm(length(tm))
% QR(1) = 0;
% for i=2:1:length(xint)
%     QR(i) = quad(f1,5e-3,xint(i));
% end



%integrate the ifcc power to find the loss
pp2 = interp1(tm,abs(PtauSUM),'linear','pp');
f2 = @(x) ppval(pp2,x);
QIFCCtot = quad(f2,5e-3,tm(length(tm)))


xint = 5e-3:5e-3:tm(length(tm))
QIFCC(1) = 0;
for i=2:1:length(xint)
    QIFCC(i) = quad(f2,5e-3,xint(i));
end


%integrate the ifcc power to find the loss
pp3 = interp1(tm,abs(PresSUM),'linear','pp');
f3 = @(x) ppval(pp3,x);
QREStot = quad(f3,5e-3,tm(length(tm)))


xint = 5e-3:5e-3:tm(length(tm))
QRES(1) = 0;
for i=2:1:length(xint)
    QRES(i) = quad(f3,5e-3,xint(i));
end


QREStot = QREStot*Lc*sym;
QRES = QRES*Lc*sym;

QIFCCtot = QIFCCtot*Li*sym;
QIFCC = QIFCC*Li*sym;


% dlmwrite('ANSYS_I_6B.csv',[tm I0])
% dlmwrite('ANSYS_ind_6B.csv',[tm LL])
% dlmwrite('ANSYS_phi_6B.csv',[tm flux])
% dlmwrite('ANSYS_Edump_6B.csv',[xint' QR'])
dlmwrite('ANSYS_Eifcc_6B.csv',[xint' QIFCC'])
dlmwrite('ANSYS_tauHT0_6B.csv',[tm tauHT0])
dlmwrite('ANSYS_MxHT0_6B.csv',[tm MxHT0])
dlmwrite('ANSYS_MyHT0_6B.csv',[tm MyHT0])
dlmwrite('ANSYS_QifccHT0_6B.csv',[tm PtauHT0])
dlmwrite('ANSYS_QifccTOT_6B.csv',[tm PtauSUM])
dlmwrite('ANSYS_Tave_HT0_6B.csv',[tmT TaveHT0])
dlmwrite('ANSYS_Thot_6B.csv',[tmTH Thot])
dlmwrite('ANSYS_Ecoil_6B.csv',[xint' QRES'])
dlmwrite('ANSYS_QresHT0_6B.csv',[tm PresHT0])
dlmwrite('ANSYS_QresTOT_6B.csv',[tm PresSUM])
dlmwrite('ANSYS_Rcoil_6B.csv',[tm Rcoil])


% datcern = importdata('D:\ansys_UDF_dist_v1\CERN_benchmarking_v11_clean\from_Lorenzo_10_8_2018\6B\all_dat_6B.txt');
datcern = importdata('..\CERN_res\6B\comsol_res_6B.txt');
tmC = datcern(:,1);
tmC = tmC + 5e-3;   %needs offset to match

tauHT0_cern = datcern(:,2);
MxHT0_cern  = datcern(:,3);
MyHT0_cern  = datcern(:,4);
PtauHT0_cern  = datcern(:,5);
PtauSUM_cern  = datcern(:,6);
PresHT0_cern  = datcern(:,7);
PresSUM_cern = datcern(:,8);
TaveHT0_cern  = datcern(:,9);
Thot_cern = datcern(:,10);
Rcoil_cern  = datcern(:,11);
Eifcc_cern  = datcern(:,12);
Ecoil_cern  = datcern(:,13);

%fix lorenzo based on email 11/29
% Dear Lucas,
% My bad. I multiplied Rcoil times 9.2[m] (magnetic length), instead of 10.11[m] (coil length).
% I fixed the calculation, and the agreement is now nearly perfect J
Rcoil_cern = Rcoil_cern*10.11/9.2;

% I0M_cern  = importdata('..\CERN_res\5D\Iext_5D.csv');
% LL_cern  = importdata('..\CERN_res\5D\L_5D.csv');
% flux_cern  = importdata('..\CERN_res\5D\phi_5D.csv');
% Edump_cern  = importdata('..\CERN_res\5D\Edump_5D.csv');
% PresHT0_cern  = importdata('..\CERN_res\5D\Q_res_tot_HT0_5D.csv');
% PresSUM_cern  = importdata('..\CERN_res\5D\Q_res_tot_5D.csv');
% Rcoil_cern  = importdata('..\CERN_res\5D\R_tot_5D.csv');
% tauHT0_cern  = importdata('..\CERN_res\5D\tau_avg_HT0_5D.csv');
% MxHT0_cern  = importdata('..\CERN_res\5D\M_ifccx_avg_HT0_5D.csv');
% MyHT0_cern  = importdata('..\CERN_res\5D\M_ifccy_avg_HT0_5D.csv');
% PtauHT0_cern  = importdata('..\CERN_res\5D\Q_ifcc_tot_HT0_5D.csv');
% PtauSUM_cern  = importdata('..\CERN_res\5D\Q_ifcc_tot_5D.csv');
% TaveHT0_cern  = importdata('..\CERN_res\5D\T_avg_HT0_5D.csv');



%check the percentage missing due to voltage supply
%integrate the dump power to find the loss
% pp2 = interp1(tm,abs(PS),'linear','pp');
% f2 = @(x) ppval(pp2,x);
% QStot = quad(f2,5e-3,1)  %about 405 J - negligible



fnum=fnum+1;
h(fnum)=figure;
scatter(tm,IFCU,'filled')   %plot as negative to match ..
hold on 
% scatter(flux_cern(:,1),flux_cern(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('IFCU','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['IFCU_vs_t_',tag],'-r300')
hold off



fnum=fnum+1;
h(fnum)=figure;
scatter(tm,tauHT0,'filled')   %plot as negative to match ..
hold on 
scatter(tmC,tauHT0_cern/2.0,10,'filled')
% scatter(sdatPhi(:,1),sdatPhi(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('ave tau ifcc HT0 (s)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['tau_vs_t_',tag],'-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(tm,MxHT0,'filled')   %plot as negative to match ..
hold on 
scatter(tmC,MxHT0_cern,10,'filled')
% scatter(tm3,MxHT03,'filled')   %plot as negative to match ..
% scatter(sdatPhi(:,1),sdatPhi(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('ave Mx iffc HT0 (A/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['MxHT0_vs_t_',tag],'-r300')
hold off

fnum=fnum+1;
h(fnum)=figure;
scatter(tm,MyHT0,'filled')   %plot as negative to match ..
hold on 
scatter(tmC,MyHT0_cern,10,'filled')
% scatter(tm3,MyHT03,'filled')   %plot as negative to match ..
% scatter(sdatPhi(:,1),sdatPhi(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('ave My iffc HT0 (A/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['MyHT0_vs_t_',tag],'-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(tm,PtauHT0,'filled')   %plot as negative to match ..
hold on 
scatter(tmC,PtauHT0_cern,10,'filled')
% scatter(tm3,PtauHT03,'filled')   %plot as negative to match ..
% % scatter(sdatPhi(:,1),sdatPhi(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Q iffc HT0 (W)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['QifccHT0_vs_t_',tag],'-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(tm,PtauSUM,'filled')   %plot as negative to match ..
hold on 
scatter(tmC,PtauSUM_cern,10,'filled')
% scatter(tm3,PtauSUM3,'filled')   %plot as negative to match ..
% scatter(sdatPhi(:,1),sdatPhi(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Q iffc total - quadrant (W)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['QifccTOT_vs_t_',tag],'-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(tmT,TaveHT0,'filled')   %plot as negative to match ..
hold on 
scatter(tmC,TaveHT0_cern,10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('temp (K)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['Tave_HT0_vs_t_',tag],'-r300')
hold off

fnum=fnum+1;
h(fnum)=figure;
scatter(tmTH,Thot,'filled')   %plot as negative to match ..
hold on 
scatter(tmC,Thot_cern,10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('hotspot temp (K)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['Thot_vs_t_',tag],'-r300')
hold off



fnum=fnum+1;
h(fnum)=figure;
scatter(tm,PresHT0,'filled')   %plot as negative to match ..
hold on 
scatter(tmC,PresHT0_cern,10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Q res HT0 - quad (W/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['QresHT0_vs_t_',tag],'-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(tm,PresSUM,'filled')   %plot as negative to match ..
hold on 
scatter(tmC,PresSUM_cern,10,'filled')
% scatter(sdatPhi(:,1),sdatPhi(:,2),10,'filled')
% plot(tm,Rcoil.*I0.*I0/(Lc*4))
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Q res total - quad (W/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['QresTOT_vs_t_',tag],'-r300')
hold off



fnum=fnum+1;
h(fnum)=figure;
scatter(tm,Rcoil,'filled')   %plot as negative to match ..
hold on 
scatter(tmC,Rcoil_cern,10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('R coil (ohm)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['Rcoil_vs_t_',tag],'-r300')
hold off

% fnum=fnum+1;
% h(fnum)=figure;
% scatter(xint,QR,'filled')   
% hold on 
% % scatter(Edump_cern(:,1),Edump_cern(:,2),10,'filled')
% scatter(xint,QIFCC,'filled') 
% scatter(xint,QRES,'filled')
% % scatter(sdatPhi(:,1),sdatPhi(:,2),10,'filled')
% plot([0 tm(length(tm))],[dE dE],'linewidth',2)
% % plot([0 1],[dE*Li/Lc dE*Li/Lc],'linewidth',2)
% box on
% hold on
% grid on
% title(tag)
% xlabel('tm (s)','FontSize',18)
% ylabel('Energy (J)','FontSize',18)
% set(gca,'FontSize',16,'linewidth',2)
% set(h(fnum),'Position', [200 200 850 600])
% legend('E dump','E dump cern','E IFCC (J)','E Coil (J)','dE 5ms to 500ms','Location','NorthEastOutside')
% set(gcf,'PaperPositionMode','auto')
% print(h(fnum),'-djpeg',['E_dump_vs_t_',tag],'-r300')
% hold off



return 

