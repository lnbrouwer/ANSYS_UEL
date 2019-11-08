clearvars 
close all
fnum = 0;


tag = 'Test0A';
dat = importdata('./user102/Mag_v_t.txt');


sym = 4;  
I0 = dat(:,1);
flux = dat(:,2)*sym;
LL = dat(:,3)*sym;
BaveHT0 = dat(:,4);
tm = dat(:,5);
Mxave = dat(:,6);
Myave = dat(:,7);
Jtauave = dat(:,8);
JtauaveALL = dat(:,9)*sym;

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


sBave = importdata('..\CERN_res\0A\B_avg_HT0_0A.csv');
sIext= importdata('..\CERN_res\0A\Iext_0A.csv');
sdatLd = importdata('..\CERN_res\0A\Ld_0A.csv');
sdatPhi = importdata('..\CERN_res\0A\phi_0A.csv');



fnum=fnum+1;
h(fnum)=figure;
scatter(tm,I0,'filled')
hold on 
scatter(sIext(:,1),sIext(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Iext (A)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['Iext_vs_t_',tag],'-r300')
hold off

fnum=fnum+1;
h(fnum)=figure;
scatter(tm,BaveHT0,'filled')
hold on 
scatter(sBave(:,1),sBave(:,2),10,'filled')
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
print(h(fnum),'-djpeg',['BaveHT0_vs_t_',tag],'-r300')
hold off



fnum=fnum+1;
h(fnum)=figure;
scatter(tm,flux,'filled')
hold on 
scatter(sdatPhi(:,1),sdatPhi(:,2),10,'filled')
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
scatter(tm,dLL,'filled')
hold on 
scatter(sdatLd(:,1),sdatLd(:,2),10,'filled')
box on
hold on
grid on
title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('Ld (H/m)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg',['Ld_vs_t_',tag],'-r300')
hold off

dlmwrite('ANSYS_B_avg_HT0_0A.csv',[tm BaveHT0])
dlmwrite('ANSYS_Iext_0A.csv',[tm I0])
dlmwrite('ANSYS_Ld_0A.csv',[tm dLL])
dlmwrite('ANSYS_phi_0A.csv',[tm flux])

return 






