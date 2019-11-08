clearvars 
close all
fnum = 0;


tag = 'Test3Bc';

dat2 = importdata('./user102/Therm_v_t.txt');
tmT = dat2(:,1);
TaveHT0 = dat2(:,2);
BaveHT0 = dat2(:,3);
CvCuHT0 = dat2(:,4);
CvNb3SnHT0 = dat2(:,5);
CvG10HT0 = dat2(:,6);
CvmixHT0 = dat2(:,7);
kCuHT0 = dat2(:,8);
kmixHT0 = dat2(:,9);



dlmwrite('ANSYS_T_avg_HT0_3Bc.csv',[tmT TaveHT0])
dlmwrite('ANSYS_B_avg_HT0_3Bc.csv',[tmT BaveHT0])
dlmwrite('ANSYS_CvCu_avg_HT0_3Bc.csv',[tmT CvCuHT0])
dlmwrite('ANSYS_Nb3Sn_avg_HT0_3Bc.csv',[tmT CvNb3SnHT0])
dlmwrite('ANSYS_G10_avg_HT0_3Bc.csv',[tmT CvG10HT0])
dlmwrite('ANSYS_Cvmix_avg_HT0_3Bc.csv',[tmT CvmixHT0])
dlmwrite('ANSYS_kmix_avg_HT0_3Bc.csv',[tmT kmixHT0])
dlmwrite('ANSYS_kCu_avg_HT0_3Bc.csv',[tmT kCuHT0])



TaveHT0_cern = importdata('..\..\CERN_res\3Bc\1e6\T_avg_HT0_3Bc.csv');
CvCuHT0_cern = importdata('..\..\CERN_res\3Bc\1e6\CvCu_avg_HT0_3Bc.csv');
CvNb3SnHT0_cern = importdata('..\..\CERN_res\3Bc\1e6\CvSc_avg_HT0_3Bc.csv');
CvG10HT0_cern = importdata('..\..\CERN_res\3Bc\1e6\CvG10_avg_HT0_3Bc.csv');
CvmixHT0_cern = importdata('..\..\CERN_res\3Bc\1e6\Cv_avg_HT0_3Bc.csv');
kmixHT0_cern = importdata('..\..\CERN_res\3Bc\1e6\k_avg_HT0_3Bc.csv');
kCuHT0_cern = importdata('..\..\CERN_res\3Bc\1e6\kCu_avg_HT0_3Bc.csv');




fnum=fnum+1;
h(fnum)=figure;
scatter(tmT,TaveHT0,'filled')
hold on 
scatter(TaveHT0_cern(:,1),TaveHT0_cern(:,2),10,'filled')
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


fnum=fnum+1;
h(fnum)=figure;
scatter(tmT,BaveHT0,'filled')
hold on 
% scatter(TaveHT0_cern(:,1),TaveHT0_cern(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('ave B HT0 (T)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['B_HT0_vs_t_',tag],'-r300')
hold off

fnum=fnum+1;
h(fnum)=figure;
scatter(tmT,CvCuHT0,'filled')
hold on 
scatter(CvCuHT0_cern(:,1),CvCuHT0_cern(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('ave cvCu HT0 (J/kg*m^3)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['cvCu_HT0_vs_t_',tag],'-r300')
hold off



fnum=fnum+1;
h(fnum)=figure;
scatter(tmT,CvNb3SnHT0,'filled')
hold on 
scatter(CvNb3SnHT0_cern(:,1),CvNb3SnHT0_cern(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('ave cvNb3Sn HT0 (J/kg*m^3)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['cvNb3Sn_HT0_vs_t_',tag],'-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(tmT,CvG10HT0,'filled')
hold on 
scatter(CvG10HT0_cern(:,1),CvG10HT0_cern(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('ave cvG10 HT0 (J/kg*m^3)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['cvG10_HT0_vs_t_',tag],'-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(tmT,CvmixHT0,'filled')
hold on 
scatter(CvmixHT0_cern(:,1),CvmixHT0_cern(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('ave cv mix HT0 (J/kg*m^3)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['cvmix_HT0_vs_t_',tag],'-r300')
hold off




fnum=fnum+1;
h(fnum)=figure;
scatter(TaveHT0,CvCuHT0,'filled')
hold on 
scatter(TaveHT0_cern(:,2),CvCuHT0_cern(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('temp (K)','FontSize',18)
ylabel('ave cvCu HT0 (J/kg*m^3)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['cvCu_HT0_vs_temp_',tag],'-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(TaveHT0,CvNb3SnHT0,'filled')
hold on 
scatter(TaveHT0_cern(:,2),CvNb3SnHT0_cern(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('temp (K)','FontSize',18)
ylabel('ave cvNb3Sn HT0 (J/kg*m^3)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['cvNb3Sn_HT0_vs_temp_',tag],'-r300')
hold off

fnum=fnum+1;
h(fnum)=figure;
scatter(TaveHT0,CvG10HT0,'filled')
hold on 
scatter(TaveHT0_cern(:,2),CvG10HT0_cern(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('temp (K)','FontSize',18)
ylabel('ave cvG10 HT0 (J/kg*m^3)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['cvG10_HT0_vs_temp_',tag],'-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(TaveHT0,CvmixHT0,'filled')
hold on 
scatter(TaveHT0_cern(:,2),CvmixHT0_cern(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('temp (K)','FontSize',18)
ylabel('ave cv mix HT0 (J/kg*m^3)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['cvmix_HT0_vs_temp_',tag],'-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(TaveHT0,kmixHT0,'filled')
hold on 
scatter(TaveHT0_cern(:,2),kmixHT0_cern(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('temp (K)','FontSize',18)
ylabel('ave k mix HT0 (W/mK)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['kmix_HT0_vs_temp_',tag],'-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(TaveHT0,kCuHT0,'filled')
hold on 
scatter(TaveHT0_cern(:,2),kCuHT0_cern(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('temp (K)','FontSize',18)
ylabel('ave k mix HT0 (W/mK)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['kCu_HT0_vs_temp_',tag],'-r300')
hold off




return


fnum=fnum+1;
h(fnum)=figure;
scatter(tm,flux,'filled')
hold on 
scatter(sdatPhi(:,1),sdatPhi(:,2),10,'filled')
box on
hold on
grid on
xlabel('tm (s)','FontSize',18)
ylabel('Linked Flux (Wb/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['flux_vs_t',tag],'-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(tm,dLL,'filled')
hold on 
scatter(sdatLd(:,1),sdatLd(:,2),10,'filled')
box on
hold on
grid on
xlabel('tm (s)','FontSize',18)
ylabel('Ld (H/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['Ld_vs_t',tag],'-r300')
hold off

dlmwrite('ANSYS_B_avg_HT0_0B.csv',[tm BaveHT0])
dlmwrite('ANSYS_Iext_0B.csv',[tm I0])
dlmwrite('ANSYS_Ld_0B.csv',[tm dLL])
dlmwrite('ANSYS_phi_0B.csv',[tm flux])

return 

dlmwrite('test_0A_user102.txt',[tm(1:length(tm)) I0(1:length(I0)) flux(1:length(flux)) LL(1:length(LL)) dLL'])
    

return 






