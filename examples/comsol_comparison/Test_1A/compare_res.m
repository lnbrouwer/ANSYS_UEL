clearvars 
close all
fnum = 0;


tag = 'Test1A';
dat = importdata('./user102/Mag_v_t.txt');


sym = 4;  
I0 = dat(:,1);
flux = dat(:,2)*sym;
LL = dat(:,3)*sym;
BaveHT0 = dat(:,4);
tm = dat(:,5);
Mxave = dat(:,6);
Myave = dat(:,7);
JtauHT0 = dat(:,8);
JtauALL = dat(:,9)*sym;

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


sMxHT0 = importdata('..\CERN_res\1A\M_ifccx_avg_HT0_1A.csv');
sMyHT0 = importdata('..\CERN_res\1A\M_ifccy_avg_HT0_1A.csv');
sQtot = importdata('..\CERN_res\1A\Q_ifcc_tot_1A.csv');
sQHT0 = importdata('..\CERN_res\1A\Q_ifcc_tot_HT0_1A.csv');


fnum=fnum+1;
h(fnum)=figure;
scatter(tm,Mxave,'filled')
hold on 
scatter(sMxHT0(:,1),sMxHT0(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Mxave HT0 (A/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['Mx_ave_HT0_vs_t_',tag],'-r300')
hold off

fnum=fnum+1;
h(fnum)=figure;
scatter(tm,Myave,'filled')
hold on 
scatter(sMyHT0(:,1),sMyHT0(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Myave HT0 (A/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['My_ave_HT0_vs_t_',tag],'-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(tm,JtauALL,'filled')
hold on 
scatter(sQtot(:,1),sQtot(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Qtot IFCC (W/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['Qtot_IFCC_vs_t_',tag],'-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(tm,JtauHT0,'filled')
hold on 
scatter(sQHT0(:,1),sQHT0(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('QHT0 IFCC (W/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['QHT0_IFCC_vs_t_',tag],'-r300')
hold off

% fnum=fnum+1;
% h(fnum)=figure;
% scatter(tm,I0,'filled')
% hold on 
% box on
% hold on
% grid on
% title(tag)
% xlabel('tm (s)','FontSize',18)
% ylabel('I (A) ','FontSize',18)
% set(gca,'FontSize',16,'linewidth',2)
% set(h(fnum),'Position', [200 200 850 600])
% legend('ANSYS','COMSOL','Location','NorthEastOutside')
% set(gcf,'PaperPositionMode','auto')
% print(h(fnum),'-djpeg',['I0_vs_t_',tag],'-r300')
% hold off


dlmwrite('ANSYS_M_ifccx_avg_HT0_1A.csv',[tm Mxave])
dlmwrite('ANSYS_M_ifccy_avg_HT0_1A.csv',[tm Myave])
dlmwrite('ANSYS_Q_ifcc_tot_1A.csv',[tm JtauALL])
dlmwrite('ANSYS_Q_ifcc_tot_HT0_1A.csv',[tm JtauHT0])

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






