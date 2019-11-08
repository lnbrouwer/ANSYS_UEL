clearvars 
close all
fnum = 0;

%convert circu res
run .\user102\write_trans_elem_res.m
datCirc = importdata('./user102/circu_res.txt');    %dlmwrite('circu_res.txt',[tm I1 V1 VR VS]);
I0 = datCirc(:,2);
VC = datCirc(:,3);
VR = datCirc(:,4);
VS = datCirc(:,5);
PC = datCirc(:,6);
PR = datCirc(:,7);
PS = datCirc(:,8);
tag = 'Test4C';

dat = importdata('./user102/Mag_v_t.txt');
% dat2 = importdata('./user102/Therm_v_t.txt');
% tmT = dat2(:,1);
% TaveHT0 = dat2(:,2);


% 
sym = 4;  
Li = 9.20
Lc = 10.11

I0M = abs(I0);
flux = dat(:,2)*sym;
LL = flux./I0;    %now divide by current is seperate (circuit)
EE = dat(:,3);
% Ei = EE(1)*sym*Li;
% % Ei = 22649*sym*Li;  %senergy for all but conductor using senergy macro
% % Ei = 26542.8*sym*Li;  %senergy for all with yoke done using macro
% Ef = EE(length(EE))*sym*Li;
Ei = (EE(1)+dat(1,4))*sym*Li;
Ef = (EE(length(EE))+dat(1,5))*sym*Li;
dE = Ei-Ef;
Emiss = Ef/Ei;

% tauHT0 = dat(:,6);
% MxHT0 = dat(:,7);
% MyHT0 = dat(:,8);
% PtauHT0 = dat(:,9);
% PtauSUM = dat(:,10);
PresHT0 = dat(:,11);
PresSUM = dat(:,12);
Rcoil = dat(:,13)*sym;

%calculate Diff L
for i=1:1:length(I0)
    if i==1
        dphi(i) = flux(i)-0;
        di(i) = I0(i)-0;
    %     dt(i-1) = tm(i)-tm(i-1);
        dLL(i) = dphi(i)/di(i);          
    else
        dphi(i) = flux(i)-flux(i-1);
        di(i) = I0(i)-I0(i-1);
    %     dt(i-1) = tm(i)-tm(i-1);
        dLL(i) = dphi(i)/di(i);        
    end
end
dLL = dLL';


%integrate the dump power to find the loss
pp = interp1(tm,abs(PR),'linear','pp');
f1 = @(x) ppval(pp,x);
QRtot = quad(f1,5e-3,tm(length(tm)))


xint = 5e-3:5e-3:tm(length(tm))
QR(1) = 0;
for i=2:1:length(xint)
    QR(i) = quad(f1,5e-3,xint(i));
end



% %integrate the ifcc power to find the loss
% pp2 = interp1(tm,abs(PtauSUM),'linear','pp');
% f2 = @(x) ppval(pp2,x);
% QIFCCtot = quad(f2,5e-3,tm(length(tm)))
% 
% 
% xint = 5e-3:5e-3:tm(length(tm))
% QIFCC(1) = 0;
% for i=2:1:length(xint)
%     QIFCC(i) = quad(f2,5e-3,xint(i));
% end

% QIFCCtot = QIFCCtot*Li*sym;
% QIFCC = QIFCC*Li*sym;

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




dlmwrite('ANSYS_I_4C.csv',[tm I0])
dlmwrite('ANSYS_ind_4C.csv',[tm LL])
dlmwrite('ANSYS_phi_4C.csv',[tm flux])
dlmwrite('ANSYS_Edump_4C.csv',[xint' QR'])
dlmwrite('ANSYS_Ecoil_4C.csv',[xint' QRES'])
dlmwrite('ANSYS_QresHT0_4C.csv',[tm PresHT0])
dlmwrite('ANSYS_QresTOT_4C.csv',[tm PresSUM])
dlmwrite('ANSYS_Rcoil_4C.csv',[tm Rcoil])


I0M_cern  = importdata('..\CERN_res\4C\Iext_4C.csv');
LL_cern  = importdata('..\CERN_res\4C\L_4C.csv');
flux_cern  = importdata('..\CERN_res\4C\phi_4C.csv');
Edump_cern  = importdata('..\CERN_res\4C\Edump_4C.csv');
PresHT0_cern  = importdata('..\CERN_res\4C\Q_res_tot_HT0_4C.csv');
PresSUM_cern  = importdata('..\CERN_res\4C\Q_res_tot_4C.csv');
Rcoil_cern  = importdata('..\CERN_res\4C\R_tot_4C.csv');



%check the percentage missing due to voltage supply
%integrate the dump power to find the loss
% pp2 = interp1(tm,abs(PS),'linear','pp');
% f2 = @(x) ppval(pp2,x);
% QStot = quad(f2,5e-3,1)  %about 405 J - negligible



fnum=fnum+1;
h(fnum)=figure;
scatter(tm,I0M,'filled')
hold on 
scatter(I0M_cern(:,1),I0M_cern(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('I (A)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['I_vs_t_',tag],'-r300')
hold off

fnum=fnum+1;
h(fnum)=figure;
scatter(tm,-flux,'filled')   %plot as negative to match ..
hold on 
scatter(flux_cern(:,1),flux_cern(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Linked Flux (Wb/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['flux_vs_t_',tag],'-r300')
hold off



fnum=fnum+1;
h(fnum)=figure;
scatter(tm,PresHT0,'filled')   %plot as negative to match ..
hold on 
scatter(PresHT0_cern(:,1),PresHT0_cern(:,2),10,'filled')
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
% scatter(sdatPhi(:,1),sdatPhi(:,2),10,'filled')
% plot(tm,Rcoil.*I0.*I0/(Lc*4))
scatter(PresSUM_cern(:,1),PresSUM_cern(:,2)/sym,10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Q res total - quad (W/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
% legend('ANSYS','check with I*I*rcoil/Lc*4','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['QresTOT_vs_t_',tag],'-r300')
hold off



fnum=fnum+1;
h(fnum)=figure;
scatter(tm,Rcoil,'filled')   %plot as negative to match ..
hold on 
scatter(Rcoil_cern(:,1),Rcoil_cern(:,2),10,'filled')
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
% scatter(tm,PR,'filled')   %plot as negative to match ..
% hold on 
% scatter(tm,P1,'filled')   %plot as negative to match ..
% scatter(tm,PS,'filled')   %plot as negative to match ..
% % scatter(sdatPhi(:,1),sdatPhi(:,2),10,'filled')
% plot(tm,f1(tm))
% % plot(tm,(25e-3)*I0.*I0)
% box on
% hold on
% grid on
% title(tag)
% xlabel('tm (s)','FontSize',18)
% ylabel('Q Dump Res. (J)','FontSize',18)
% set(gca,'FontSize',16,'linewidth',2)
% set(h(fnum),'Position', [200 200 850 600])
% legend('ANSYS','COMSOL','Location','NorthEastOutside')
% set(gcf,'PaperPositionMode','auto')
% print(h(fnum),'-djpeg',['flux_vs_t_',tag],'-r300')
% hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(xint,QR,'filled')   
hold on 
% scatter(xint,QIFCC,'filled') 
scatter(Edump_cern(:,1),Edump_cern(:,2),10,'filled')
scatter(xint,QRES,'filled')
% scatter(sdatPhi(:,1),sdatPhi(:,2),10,'filled')
plot([0 tm(length(tm))],[dE dE],'linewidth',2)
% plot([0 1],[dE*Li/Lc dE*Li/Lc],'linewidth',2)
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Energy (J)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('E dump','E dump Cern','E RES (J)','dE 5ms to 500ms','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['E_dump_vs_t_',tag],'-r300')
hold off




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compare with other test vs. I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dat = importdata('../Test_0B/user102/Mag_v_t.txt');
I2 = dat(:,1);
flux2 = dat(:,2)*sym;
LL2 = dat(:,3)*sym;


fnum=fnum+1;
h(fnum)=figure;
scatter(I0M,LL,'filled')
hold on 
scatter(I0M_cern(:,2),LL_cern(:,2),10,'filled')
plot(I2,LL2,'linewidth',2)
box on
hold on
grid on
title(tag)
xlabel('I0 (A)','FontSize',18)
ylabel('L (H/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('4C ANSYS','4C CERN','0B','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','LL_vs_I_compare','-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(I0M,-flux,'filled')
hold on 
scatter(I0M_cern(:,2),flux_cern(:,2),10,'filled')
plot(I2,flux2,'linewidth',2)
% scatter(sdatLd(:,1),sdatLd(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('I0 (A)','FontSize',18)
ylabel('Linked Flux (Wb/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('4C ANSYS','4C CERN','0B','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','flux_vs_I_compare','-r300')
hold off

