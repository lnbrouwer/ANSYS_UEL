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
tag = 'Test5D';

dat = importdata('./user102/Mag_v_t.txt');
dat2 = importdata('./user102/Therm_v_t.txt');
tmT = dat2(:,1);
TaveHT0 = dat2(:,2);

dat3 = importdata('./user102/Thot_v_t.txt');
tmTH = dat3(:,1);
Thot = dat3(:,2);
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


tauHT0 = dat(:,6);
MxHT0 = dat(:,7);
MyHT0 = dat(:,8);
PtauHT0 = dat(:,9);
PtauSUM = dat(:,10);

PresHT0 = dat(:,11);
PresSUM = dat(:,12);
Rcoil = dat(:,13)*sym;
IFCU = dat(:,14);



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


dlmwrite('ANSYS_I_5D.csv',[tm I0])
dlmwrite('ANSYS_ind_5D.csv',[tm LL])
dlmwrite('ANSYS_phi_5D.csv',[tm flux])
dlmwrite('ANSYS_Edump_5D.csv',[xint' QR'])
dlmwrite('ANSYS_Eifcc_5D.csv',[xint' QIFCC'])
dlmwrite('ANSYS_tauHT0_5D.csv',[tm tauHT0])
dlmwrite('ANSYS_MxHT0_5D.csv',[tm MxHT0])
dlmwrite('ANSYS_MyHT0_5D.csv',[tm MyHT0])
dlmwrite('ANSYS_QifccHT0_5D.csv',[tm PtauHT0])
dlmwrite('ANSYS_QifccTOT_5D.csv',[tm PtauSUM])
dlmwrite('ANSYS_Tave_HT0_5D.csv',[tmT TaveHT0])
dlmwrite('ANSYS_Thot_5D.csv',[tmTH Thot])
dlmwrite('ANSYS_Ecoil_5D.csv',[xint' QRES'])
dlmwrite('ANSYS_QresHT0_5D.csv',[tm PresHT0])
dlmwrite('ANSYS_QresTOT_5D.csv',[tm PresSUM])
dlmwrite('ANSYS_Rcoil_5D.csv',[tm Rcoil])

dlmwrite('ANSYS_Rcoil_5D.csv',[TaveHT0 Rcoil])


%now import all at once
datc = importdata('..\CERN_res\5D\comsol_res_5D.txt');
tmc = datc(:,1);
I0M_cern  = datc(:,2);
LL_cern  = datc(:,4);
flux_cern = datc(:,3);
Edump_cern  = datc(:,15);
PresHT0_cern  = datc(:,10);
PresSUM_cern  = datc(:,11);
Rcoil_cern  = datc(:,14);
tauHT0_cern  = datc(:,5);
MxHT0_cern  = datc(:,6);
MyHT0_cern  = datc(:,7);
PtauHT0_cern  = datc(:,8);
PtauSUM_cern  = datc(:,9);
TaveHT0_cern  = datc(:,12);
Thot_cern = datc(:,13);



%CERN DATA NEEDS TO BE SHIFTED BY 5MS (it starts at zero)
tmc =  tmc  + 5e-3;

%there is a duplicate row at time = 0.02
% [~,idx] = unique(PtauSUM_cern);   %which rows have a unique first value?
% PtauSUM_cern = PtauSUM_cern(idx,:)                %only use those

%integrate the ifcc power to find the loss
pp2 = interp1(tmc,abs(PtauSUM_cern),'linear','pp');
f2 = @(x) ppval(pp2,x);
QIFCCtot_CERN = quad(f2,5e-3,tmc(length(tmc)))
QIFCCtot_CERN = QIFCCtot_CERN*Li*sym;


xintC = 5e-3:5e-3:tmc(length(tmc))
QIFCC_CERN(1) = 0;
for i=2:1:length(xintC)
    QIFCC_CERN(i) = quad(f2,5e-3,xintC(i));
end
QIFCC_CERN = QIFCC_CERN*Li*sym;



%integrate the ifcc power to find the loss
pp3 = interp1(tmc,abs(PresSUM_cern),'linear','pp');
f3 = @(x) ppval(pp3,x);
QREStot_cern = quad(f3,5e-3,tmc(length(tmc)))
QRES_cern(1) = 0;
for i=2:1:length(xintC)
    QRES_cern(i) = quad(f3,5e-3,xintC(i));
end
QREStot_cern = QREStot_cern*Lc*sym;
QRES_cern = QRES_cern*Lc*sym;



%comparison to other IFCC case with no temp
dat3 = importdata('../Test_5A/user102/Mag_v_t.txt');

tm3 = dat3(:,1);
tauHT03 = dat3(:,6);
MxHT03 = dat3(:,7);
MyHT03 = dat3(:,8);
PtauHT03 = dat3(:,9);
PtauSUM3 = dat3(:,10);


%check the percentage missing due to voltage supply
%integrate the dump power to find the loss
% pp2 = interp1(tm,abs(PS),'linear','pp');
% f2 = @(x) ppval(pp2,x);
% QStot = quad(f2,5e-3,1)  %about 405 J - negligible



fnum=fnum+1;
h(fnum)=figure;
scatter(tm,I0M,'filled')
hold on 
scatter(tmc,I0M_cern,10,'filled')
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
scatter(tmc,flux_cern,10,'filled')
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
scatter(tmc,tauHT0_cern/2.0,10,'filled')
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
scatter(tmc,MxHT0_cern,10,'filled')
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
scatter(tmc,MyHT0_cern,10,'filled')
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
scatter(tmc,PtauHT0_cern,10,'filled')
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
scatter(tmc,PtauSUM_cern,10,'filled')
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
scatter(tmc,TaveHT0_cern,10,'filled')
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
scatter(tmc,Thot_cern,10,'filled')
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
scatter(tmc,PresHT0_cern,10,'filled')
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
scatter(tmc,PresSUM_cern,10,'filled')
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
scatter(tmc,Rcoil_cern,10,'filled')
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

fnum=fnum+1;
h(fnum)=figure;
scatter(xint,QR,'filled')   
hold on 
scatter(tmc,Edump_cern,10,'filled')
scatter(xint,QIFCC,'filled') 
scatter(xintC,QIFCC_CERN,10,'filled') 
scatter(xint,QRES,'filled')
scatter(xintC,QRES_cern,10,'filled') 
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
legend('E dump','E dump cern','E IFCC','E IFCC cern','E Coil','E Coil cern','dE 5ms to 500ms','Location','NorthEastOutside')
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
% scatter(I0M_cern(:,2),LL_cern(:,2),10,'filled')
plot(I2,LL2,'linewidth',2)
% scatter(sdatLd(:,1),sdatLd(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('I0 (A)','FontSize',18)
ylabel('L (H/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('5D ANSYS','5D COMSOL','0B','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','LL_vs_I_compare','-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(I0M,-flux,'filled')
hold on 
% scatter(I0M_cern(:,2),flux_cern(:,2),10,'filled')
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
legend('5D ANSYS','5D COMSOL','0B','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','flux_vs_I_compare','-r300')
hold off

%compare linked flux vs. I


return 


fnum=fnum+1;
h(fnum)=figure;
scatter(tm,VC,'filled')
hold on 
scatter(tm,VR,'filled')
scatter(tm,VS,'filled')
% scatter(sTHT0(:,1),sTHT0(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('V (volt)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('Vcoil','Vdump','Vsource','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['V_vs_t_',tag],'-r300')
hold off



fnum=fnum+1;
h(fnum)=figure;
scatter(tm,JsumHT0,'filled')
hold on 
% scatter(sQHT0(:,1),sQHT0(:,2),10,'filled')
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
% scatter(sQtot(:,1),sQtot(:,2),10,'filled')
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
% scatter(sRhoHT0(:,1),sRhoHT0(:,2),10,'filled')
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
scatter(tm,Rcoil*1e3,'filled')
hold on 
% scatter(sRhoHT0(:,1),sRhoHT0(:,2),10,'filled')
% axis([sRhoHT0(1,1) sRhoHT0(length(sRhoHT0(:,1)),1) 0.95*min(sRhoHT0(:,2)) 1.05*max(sRhoHT0(:,2))])
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Rcoil (mohm)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['Rcoil_vs_t_',tag],'-r300')
hold off




fnum=fnum+1;
h(fnum)=figure;
scatter(tmT,TaveHT0,'filled')
hold on 
% scatter(sTHT0(:,1),sTHT0(:,2),10,'filled')
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


%%%%use calculated Ld to compare to expected LR decay 




dlmwrite('ANSYS_Q_avg_HT0_3A_CC.csv',[tm JsumHT0])
dlmwrite('ANSYS_Q_tot_3A_CC.csv',[tm JALL])
dlmwrite('ANSYS_rho_avg_HT0_3A_CC.csv',[tm rhoHT0])
dlmwrite('ANSYS_T_avg_HT0_3A_CC.csv',[tm TaveHT0])

return 


%estimate coil resistance (see if it is anywhere close to dump) (should this account for #of turns somehow?)
rhoave = 1.2e-9;
sc = (0.36e-3);   
lc = 1.0;
nc = 16;
Rcoil = (rhoave*lc/sc)*1e3*sym  %in mOhm

Rcoil = (lc*(nc/sc)*(nc/sc)*rhoave*sc)*1e3*sym  %in mOhm



nc = 1:1:16
sc1 = 0.22500e-4;
sc = sc1*nc;
Rcoil = (lc*(nc./sc).*(nc./sc)*rhoave.*sc)*1e3*sym  %in mOhm (linearly proportional to the number of turns!!) (linearly inversly proportional to area of single turn)

figure
plot(nc,Rcoil)


% ti = 5.1e-3+1e-5;
% tf = 5.5e-3;
% t = ti:1e-5:tf;
% 
% figure
% 
% plot(t,1-(exp(tf-t)/exp(tf-ti)))
% 
% l = length(t)
% exp(tf-t(1))/exp(tf-ti)
% exp(tf-t(l))/exp(tf-ti)
% 
% 
% exp(tf-t(1))
% 
% plot( exp( (tf-t))))
% 
% 
% 
% plot( exp((tf-ti)) - exp( (tf-t))))
% 
% 
% x = (t-ti)/(tf-ti);
% 
% figure
% plot(t,x.^.5)
% hold on 
% plot(t,x.^1)
% plot(t,x.^2)
