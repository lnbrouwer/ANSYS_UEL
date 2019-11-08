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

tag = 'Test6B';

dlmwrite('ANSYS_I1_6B.csv',[tm I1])
dlmwrite('ANSYS_I2_6B.csv',[tm I2])
dlmwrite('ANSYS_ICLIQ_6B.csv',[tm IRC])
dlmwrite('ANSYS_V1_6B.csv',[tm V1])
dlmwrite('ANSYS_V2_6B.csv',[tm V2])
dlmwrite('ANSYS_VCAP_6B.csv',[tm VC])


datcern = importdata('..\CERN_res\6B\comsol_res_6B_circu.txt');
tmC = datcern(:,1);
tmC = tmC + 5e-3;   %needs offset to match

I1_cern = datcern(:,2);
I2_cern  = datcern(:,3);
Icliq_cern  = datcern(:,4);
V1_cern  = datcern(:,5);
V2_cern  = datcern(:,6);
Vcliq_cern  = -datcern(:,7);




fnum=fnum+1;
h(fnum)=figure;
scatter(tm,V1,5,'filled')
hold on 
scatter(tm,V2,5,'filled')
scatter(tm,VC,5,'filled')
plot(tmC,V1_cern)
plot(tmC,V2_cern)
plot(tmC,Vcliq_cern)
axis([0 .05 -400 400])
% plot(tm,VR)
% scatter(tm,VSin,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('V (volt)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('V1','V2','VCap','V1-CERN','V2-CERN','VCap-CERN','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['V_vs_t_',tag],'-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(tm,I1,5,'filled')
hold on 
scatter(tm,I2,5,'filled')
% scatter(tm,I1-I2,5,'filled')
scatter(tm,IRC,5,'filled')
plot(tmC,I1_cern)
plot(tmC,I2_cern)
plot(tmC,Icliq_cern)
axis([0 .05 -4000 16000])
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('I (A)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('I1','I2','ICap','I1-CERN','I2-CERN','ICap-CERN','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['I_vs_t_',tag],'-r300')
hold off

fnum=fnum+1;
h(fnum)=figure;
scatter(tm,VC,'filled')
hold on 
scatter(tm,VRC,'filled')
scatter(tm,VSC,'filled')
% plot(tm,VR)
% scatter(tm,VSin,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('V (volt)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('VC','VRC','VSC','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['VCLIQ_vs_t_',tag],'-r300')
hold off


return 

%compare to previous calc with no coupling to external

% WM = importdata('real_mutual.txt');
NM = importdata('sm_mutual.txt');

fnum=fnum+1;
h(fnum)=figure;
scatter(tm,Iin,'filled')
hold on 
scatter(tm,Iout,'filled')
plot(NM(:,1),NM(:,2),'-k','linewidth',2)
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Iin (A)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('Iin','Iout','ILR-prev-noM12','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','Iin_vs_t_compare_old','-r300')
hold off

return 

fnum=fnum+1;
h(fnum)=figure;
scatter(tm,Vout,'filled')
hold on 
scatter(tm,VRout,'filled')
scatter(tm,VSout,5,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('V (volt)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('Vout','VRout','VSout','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['V_vs_t_out',tag],'-r300')
hold off

fnum=fnum+1;
h(fnum)=figure;
scatter(tm,Iout,5,'filled')
hold on
scatter(tm,IRout,'filled')
hold on
scatter(tm,ISout,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('I (A)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('Iout','IRout','ISout','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['I_vs_t_out',tag],'-r300')
hold off



fnum=fnum+1;
h(fnum)=figure;
scatter(tm,Iin,'filled')
hold on 
scatter(tm,Iout,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('I (A)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('Iin','Iout','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['I_vs_t_compare',tag],'-r300')
hold off


% dlmwrite('real_mutual.txt',[tm Iin Iout]);
% dlmwrite('sm_mutual.txt',[tm Iin Iout]);

WM = importdata('real_mutual.txt');
NM = importdata('sm_mutual.txt');

fnum=fnum+1;
h(fnum)=figure;
scatter(WM(:,1),WM(:,2),'filled')
hold on 
scatter(NM(:,1),NM(:,2),'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Iin (A)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('Iin','Iout','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','Iin_vs_t_compare_mutual','-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(WM(:,1),WM(:,3),'filled')
hold on 
scatter(NM(:,1),NM(:,3),'filled')
axis([5e-3 WM(length(WM(:,1)),1) 0.95*max(WM(:,3)) 1.05*max(WM(:,3))])
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Iout (A)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('Iin','Iout','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','Iout_vs_t_compare_mutual','-r300')
hold off

return 

