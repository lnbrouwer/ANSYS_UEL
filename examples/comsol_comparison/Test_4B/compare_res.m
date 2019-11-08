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
tag = 'Test4B';

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

tauHT0 = dat(:,6);
MxHT0 = dat(:,7);
MyHT0 = dat(:,8);
PtauHT0 = dat(:,9);
PtauSUM = dat(:,10);
%change to include full sym to match cern results
% PtauSUM = dat(:,10)*sym;

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


QIFCCtot = QIFCCtot*Li*sym;
QIFCC = QIFCC*Li*sym;


dlmwrite('ANSYS_I_4B.csv',[tm I0])
dlmwrite('ANSYS_ind_4B.csv',[tm LL])
dlmwrite('ANSYS_phi_4B.csv',[tm flux])
dlmwrite('ANSYS_Edump_4B.csv',[xint' QR'])
dlmwrite('ANSYS_Eifcc_4B.csv',[xint' QIFCC'])
dlmwrite('ANSYS_tauHT0_4B.csv',[tm tauHT0])
dlmwrite('ANSYS_MxHT0_4B.csv',[tm MxHT0])
dlmwrite('ANSYS_MyHT0_4B.csv',[tm MyHT0])
dlmwrite('ANSYS_QifccHT0_4B.csv',[tm PtauHT0])
dlmwrite('ANSYS_QifccTOT_4B.csv',[tm PtauSUM])


I0M_cern  = importdata('..\CERN_res\4B\Iext_4B.csv');
LL_cern  = importdata('..\CERN_res\4B\L_4B.csv');
flux_cern  = importdata('..\CERN_res\4B\phi_4B.csv');
Edump_cern  = importdata('..\CERN_res\4B\Edump_4B.csv');
tauHT0_cern  = importdata('..\CERN_res\4B\tau_avg_HT0_4B.csv');
MxHT0_cern  = importdata('..\CERN_res\4B\M_ifccx_avg_HT0_4B.csv');
MyHT0_cern  = importdata('..\CERN_res\4B\M_ifccy_avg_HT0_4B.csv');
PtauHT0_cern  = importdata('..\CERN_res\4B\Q_ifcc_tot_HT0_4B.csv');
PtauSUM_cern  = importdata('..\CERN_res\4B\Q_ifcc_tot_4B.csv');


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
scatter(tm,tauHT0,'filled')  
hold on 
scatter(tauHT0_cern(:,1),tauHT0_cern(:,2)/2.0,10,'filled')
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
scatter(MxHT0_cern(:,1),-MxHT0_cern(:,2),10,'filled')
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
scatter(tm,MyHT0,'filled')  
hold on 
scatter(MyHT0_cern(:,1),-MyHT0_cern(:,2),10,'filled')
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
scatter(PtauHT0_cern(:,1),PtauHT0_cern(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Q iffc HT0 - quadrant (W/m)','FontSize',18)
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
scatter(PtauSUM_cern(:,1),PtauSUM_cern(:,2)/sym,10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Q iffc total - quadrant (W/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['QifccTOT_vs_t_',tag],'-r300')
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
scatter(xint,QIFCC,'filled') 
scatter(Edump_cern(:,1),Edump_cern(:,2),10,'filled')
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
legend('E dump','E IFCC (J)','E dump cern','dE 5ms to 500ms','Location','NorthEastOutside')
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
legend('4B ANSYS','4B COMSOL','0B','Location','NorthEastOutside')
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
legend('4B ANSYS','4B COMSOL','0B','Location','NorthEastOutside')
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

% sQHT0 = importdata('D:\ansys_UDF_dist_v1\CERN_benchmarking_v3\COMSOL_results\3A\Q_avg_HT0_3A.csv');
% sQtot = importdata('D:\ansys_UDF_dist_v1\CERN_benchmarking_v3\COMSOL_results\3A\Q_tot_3A.csv');
% sRhoHT0 = importdata('D:\ansys_UDF_dist_v1\CERN_benchmarking_v3\COMSOL_results\3A\rho_avg_HT0_3A.csv');
% sTHT0 = importdata('D:\ansys_UDF_dist_v1\CERN_benchmarking_v3\COMSOL_results\3A\T_avg_HT0_3A.csv');


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
