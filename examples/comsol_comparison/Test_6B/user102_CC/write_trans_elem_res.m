clearvars 
close all
fnum = 0;


%read in power vs. t
% Pin = readData('power1.txt',6);
Vin = readData('voltage1.txt',6);
Iin = readData('current1.txt',6);

% set time vector - same for all
tm = Vin(:,1);
tmm = Vin(:,1)*1000;


V1 = Vin(:,2);
V2 = Vin(:,3);
VR = Vin(:,4);
VRC = Vin(:,5);
VS = Vin(:,6);
VSC = Vin(:,7);

I1 = Iin(:,2);
I2 = Iin(:,3);
IR = Iin(:,4);
IRC = Iin(:,5);
IS = Iin(:,6);
ISC = Iin(:,7);

Vin = readData('voltage2.txt',6);
Iin = readData('current2.txt',6);

VC = Vin(:,2);
IC = Vin(:,2);

dlmwrite('circu_res.txt',[tm I1 I2 IR IRC IS ISC V1 V2 VR VRC VS VSC VC IC]);

% fnum=fnum+1;
% h(fnum)=figure;
% scatter(tm,VC)
% hold on 
% scatter(tm,VR)


return 

fnum=fnum+1;
h(fnum)=figure;
scatter(tmm,abs(I1),'filled')
hold on
scatter(tmm,abs(I2),'filled')
scatter(tmm,abs(I3),'filled')
scatter(tmm,abs(I4),'filled')
scatter(tmm,abs(IR),'filled')
scatter(tmm,abs(IS),'filled')

box on
grid on
% title('Tosca In Straight Section, 5 mm ref','FontWeight','bold','FontSize',18)
xlabel('t (ms)','FontSize',18)
ylabel('I (A)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('I1','I2','I3','I4','IR','IS','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','I_vs_t','-r300')
hold off



fnum=fnum+1;
h(fnum)=figure;
scatter(tmm,V1,'filled')
hold on
scatter(tmm,V2,'filled')
scatter(tmm,V3,'filled')
scatter(tmm,V4,'filled')
scatter(tmm,VR,'filled')
scatter(tmm,VS,'filled')
box on
grid on
% title('Tosca In Straight Section, 5 mm ref','FontWeight','bold','FontSize',18)
xlabel('t (ms)','FontSize',18)
ylabel('v (volt)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('V1','V2','V3','V4','VR','VS','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','V_vs_t','-r300')
hold off




%compare to pure L/R
tmin = tm(1);
tmax = tm(length(tm));
t = tmin:1e-3:tmax

%%% calculates the stored energy based on zero freq inductance from emmanuele - per m
lmag = 0.04661;
R0 = 1.0;  %dump resistor
%L0_nofe = (8.1268e-3)*lmag;  %per m (LEDET INPUT from ROXIE)
L0_nofe = (14.3e-3)*lmag;  %per m  (LEDET INPUT FROM ANSYS IND. CALC)
I0 = IR(1);
E0 = L0_nofe*I0*I0/(2.0*10^6)


%%%%%%% Set up L/R decay comparison with roxie value
L0 = L0_nofe; 
dt = 1.0*10^-3;
t1 = 0:dt:(tmin-dt);
t2 = tmin:(1.0*10^-3):tmax;
%pre-detection
cnt = 1;
for i=1:1:length(t1)
    t(cnt) = t1(i);
    It(cnt) = I0;
    Isqt(cnt) = It(cnt).^2;
    cnt = cnt+1;
end
%extraction
for i=1:1:length(t2)
    t(cnt) = t2(i);
    It(cnt) = I0*exp(-(t2(i)-tmin)*R0/L0);
    cnt = cnt+1;
end
pp1 = interp1(t,It,'pchip','pp');
f1 = @(x) ppval(pp1,x);



%read in Nmisc for a particular element
Nin1 = readData('nmisc.txt',6);
Nin2 = readData('nmisc_2.txt',6);
Nin3 = readData('nmisc_3.txt',6);

Nin1L = readData('nmiscL.txt',6);
Nin2L = readData('nmisc_2L.txt',6);
Nin3L = readData('nmisc_3L.txt',6);
% 
% Nin1_there = readData('D:\ansys_SCU_udf\02_scu_mini\01_LR_decay\devel_v3_bat_less_msh_lmag_eng_with_eddy_keyopt3_fix\nmisc.txt',6);
% Nin2_there = readData('D:\ansys_SCU_udf\02_scu_mini\01_LR_decay\devel_v3_bat_less_msh_lmag_eng_with_eddy_keyopt3_fix\nmisc_2.txt',6);
% 


fnum=fnum+1;
h(fnum)=figure;
scatter(tmm,Nin1(:,5),5,'filled')
hold on
scatter(tmm,Nin1L(:,5),5,'filled')
% scatter(tmm,Nin1_there(:,5),5,'filled')
box on
grid on
% title('Tosca In Straight Section, 5 mm ref','FontWeight','bold','FontSize',18)
xlabel('t (ms)','FontSize',18)
ylabel('Tau (ms)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS HF','ANSYS LF','Location','NorthEastOutside')
% legend('V1O','V1I','V4I','V4O','V2I','V2O','V3O','V3I','VR','VS','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','tau_vs_t_compare_ansys','-r300')
hold off



fnum=fnum+1;
h(fnum)=figure;
scatter(tmm,Nin2(:,3),5,'filled')
hold on
scatter(tmm,Nin2L(:,3),5,'filled')
box on
grid on
% title('Tosca In Straight Section, 5 mm ref','FontWeight','bold','FontSize',18)
xlabel('t (ms)','FontSize',18)
ylabel('IFCU','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS HF','ANSYS LF','Location','NorthEastOutside')
% legend('V1O','V1I','V4I','V4O','V2I','V2O','V3O','V3I','VR','VS','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','IFCU_vs_t_compare_ansys','-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(tmm,Nin1(:,2),5,'filled')
hold on
scatter(tmm,Nin1L(:,2),5,'filled')
box on
grid on
% title('Tosca In Straight Section, 5 mm ref','FontWeight','bold','FontSize',18)
xlabel('t (ms)','FontSize',18)
ylabel('temp (K)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS HF','ANSYS LF','Location','NorthEastOutside')
% legend('V1O','V1I','V4I','V4O','V2I','V2O','V3O','V3I','VR','VS','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','temp_vs_t_compare_ansys','-r300')
hold off



fnum=fnum+1;
h(fnum)=figure;
scatter(tmm,Nin1(:,3),5,'filled')
hold on
scatter(tmm,Nin1L(:,3),5,'filled')
box on
grid on
% title('Tosca In Straight Section, 5 mm ref','FontWeight','bold','FontSize',18)
xlabel('t (ms)','FontSize',18)
ylabel('Bev (T)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS HF','ANSYS LF','Location','NorthEastOutside')
% legend('V1O','V1I','V4I','V4O','V2I','V2O','V3O','V3I','VR','VS','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','Bev_vs_t_compare_ansys','-r300')
hold off



% 
% fnum=fnum+1;
% h(fnum)=figure;
% scatter(tmm,Nin2(:,7),5,'filled')
% hold on
% % scatter(tmm,Nin2_there(:,7),5,'filled')
% box on
% grid on
% % title('Tosca In Straight Section, 5 mm ref','FontWeight','bold','FontSize',18)
% xlabel('t (ms)','FontSize',18)
% ylabel('Tau (ms)','FontSize',18)
% set(gca,'FontSize',16,'linewidth',2)
% set(h(fnum),'Position', [200 200 850 600])
% legend('ANSYS HF','ANSYS LF','ANSYS HF - there','Location','NorthEastOutside')
% % legend('V1O','V1I','V4I','V4O','V2I','V2O','V3O','V3I','VR','VS','Location','NorthEastOutside')
% set(gcf,'PaperPositionMode','auto')
% print(h(fnum),'-djpeg','tau_vs_t_compare_ansys','-r300')
% hold off


fnum=fnum+1;
h(fnum)=figure;
scatter(tmm,Nin3(:,4),5,'filled')
hold on
scatter(tmm,Nin3L(:,4),5,'filled')
% scatter(tmm,Nin3_there(:,7),5,'filled')
box on
grid on
% title('Tosca In Straight Section, 5 mm ref','FontWeight','bold','FontSize',18)
xlabel('t (ms)','FontSize',18)
ylabel('dBdt (T/s)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS HF','ANSYS LF','Location','NorthEastOutside')
% legend('V1O','V1I','V4I','V4O','V2I','V2O','V3O','V3I','VR','VS','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','dBdt_vs_t_compare_ansys','-r300')
hold off






% fnum=fnum+1;
% h(fnum)=figure;
% scatter(Nin3(:,2),Nin3(:,7),5,'filled')
% hold on
% scatter(Nin1_there(:,2),Nin2_there(:,7),5,'filled')
% box on
% grid on
% % title('Tosca In Straight Section, 5 mm ref','FontWeight','bold','FontSize',18)
% xlabel('t (ms)','FontSize',18)
% ylabel('Tau (ms)','FontSize',18)
% set(gca,'FontSize',16,'linewidth',2)
% set(h(fnum),'Position', [200 200 850 600])
% legend('ANSYS HF','ANSYS LF','ANSYS HF - there','Location','NorthEastOutside')
% % legend('V1O','V1I','V4I','V4O','V2I','V2O','V3O','V3I','VR','VS','Location','NorthEastOutside')
% set(gcf,'PaperPositionMode','auto')
% print(h(fnum),'-djpeg','temp_vs_rsvx_compare_ansys','-r300')
% hold off
% 




return 


fnum=fnum+1;
h(fnum)=figure;
plot(tmm,abs(I1))
hold on
plot(Iin2(:,1)*1000,abs(Iin2(:,2)))

box on
grid on
% title('Tosca In Straight Section, 5 mm ref','FontWeight','bold','FontSize',18)
xlabel('t (ms)','FontSize',18)
ylabel('I (A)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS with eddy','ANSYS no eddy','test','Location','NorthEastOutside')
% legend('V1O','V1I','V4I','V4O','V2I','V2O','V3O','V3I','VR','VS','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','I_vs_t_compare_ansys','-r300')
hold off






return 


Pin2 = readData4('power2.txt',6);
Pin2 = Pin2(:,2:5);
Pin = [Pin Pin2];

P1O = Pin(:,2);
P1I = Pin(:,3);
P4I = Pin(:,4);
P4O = Pin(:,5);
P2I = Pin(:,6);
P2O = Pin(:,7);
P3O = Pin(:,8);
P3I = Pin(:,9);
PR = Pin(:,10);
PS = Pin(:,11);

%read in voltage vs. t
Vin = readData('voltage1.txt',6);
Vin2 = readData4('voltage2.txt',6);
Vin2 = Vin2(:,2:5);
Vin = [Vin Vin2];

V1O = Vin(:,2);
V1I = Vin(:,3);
V4I = Vin(:,4);
V4O = Vin(:,5);
V2I = Vin(:,6);
V2O = Vin(:,7);
V3O = Vin(:,8);
V3I = Vin(:,9);
VR = Vin(:,10);
VS = Vin(:,11);

%read in current vs. t
Iin = readData('current1.txt',6);
Iin2 = readData4('current2.txt',6);
Iin2 = Iin2(:,2:5);
Iin = [Iin Iin2];

I1O = Iin(:,2);
I1I = Iin(:,3);
I4I = Iin(:,4);
I4O = Iin(:,5);
I2I = Iin(:,6);
I2O = Iin(:,7);
I3O = Iin(:,8);
I3I = Iin(:,9);
IR = Iin(:,10);
IS = Iin(:,11);


% set time vector - same for all
tm = Vin(:,1);
tmm = Vin(:,1)*1000;

fnum=fnum+1;
h(fnum)=figure;
scatter(tmm,V1O,'filled')
hold on
scatter(tmm,V1I,'filled')
scatter(tmm,V4I,'filled')
scatter(tmm,V4O,'filled')
scatter(tmm,V2I,'filled')
scatter(tmm,V2O,'filled')
scatter(tmm,V3O,'filled')
scatter(tmm,V3I,'filled')
scatter(tmm,VR,'filled')
scatter(tmm,VS,'filled')
% plot(time_vector*1000,U_CoilSections(:,1),'-r','linewidth',0.5)
% plot(time_vector*1000,U_CoilSections(:,2),'-b','linewidth',0.5)
% plot(time_vector*1000,U_CoilSections(:,8),'-m','linewidth',0.5)
% plot(time_vector*1000,U_CoilSections(:,7),'-c','linewidth',0.5)
% plot(time_vector*1000,U_CoilSections(:,4),'-k','linewidth',0.5)
% plot(time_vector*1000,U_CoilSections(:,3),'-g','linewidth',0.5)
% plot(time_vector*1000,U_CoilSections(:,5),'-b','linewidth',0.5)
% plot(time_vector*1000,U_CoilSections(:,6),'-r','linewidth',0.5)
box on
grid on
% title('Tosca In Straight Section, 5 mm ref','FontWeight','bold','FontSize',18)
xlabel('t (ms)','FontSize',18)
ylabel('V (volt)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('V1O','V1I','V4I','V4O','V2I','V2O','V3O','V3I','VR','VS','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','V_vs_t','-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
plot(tmm,I1O,'linewidth',3)
hold on
plot(tmm,I1I,'linewidth',3)
plot(tmm,I4I,'linewidth',3)
plot(tmm,I4O,'linewidth',3)
plot(tmm,I2I,'linewidth',3)
plot(tmm,I2O,'linewidth',3)
plot(tmm,I3O,'linewidth',3)
plot(tmm,I3I,'linewidth',3)
plot(tmm,IR,'linewidth',3)
% axis([0 600 0  18e3])
box on
grid on
xlabel('t (ms)','FontSize',18)
ylabel('I (A)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('I1O','I1I','I4I','I4O','I2I','I2O','I3O','I3I','IR','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','I_vs_t','-r300')
hold off


fnum=fnum+1;
h(fnum)=figure;
plot(tmm,P1O,'linewidth',3)
hold on
plot(tmm,P1I,'linewidth',3)
plot(tmm,P4I,'linewidth',3)
plot(tmm,P4O,'linewidth',3)
plot(tmm,P2I,'linewidth',3)
plot(tmm,P2O,'linewidth',3)
plot(tmm,P3O,'linewidth',3)
plot(tmm,P3I,'linewidth',3)
plot(tmm,PR,'linewidth',3)
plot(tmm,PS,'linewidth',3)
box on
grid on
xlabel('t (ms)','FontSize',18)
ylabel('P (W)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('P1O','P1I','P4I','P4O','P2I','P2O','P3O','P3I','PR','PS','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','P_vs_t','-r300')
hold off


%compare to pure L/R
tmin = tm(1);
tmax = tm(length(tm));
t = tmin:1e-3:tmax

%%% calculates the stored energy based on zero freq inductance from emmanuele - per m
lmag = 1.192;
R0 = 30e-3;  %dump resistor
%L0_nofe = (8.1268e-3)*lmag;  %per m (LEDET INPUT from ROXIE)
L0_nofe = (8.1008e-3)*lmag;  %per m  (LEDET INPUT FROM ANSYS IND. CALC)
I0 = IR(1);
E0 = L0_nofe*I0*I0/(2.0*10^6)


%%%%%%% Set up L/R decay comparison with roxie value
L0 = L0_nofe; 
dt = 1.0*10^-3;
t1 = 0:dt:(tmin-dt);
t2 = tmin:(1.0*10^-3):tmax;
%pre-detection
cnt = 1;
for i=1:1:length(t1)
    t(cnt) = t1(i);
    It(cnt) = I0;
    Isqt(cnt) = It(cnt).^2;
    cnt = cnt+1;
end
%extraction
for i=1:1:length(t2)
    t(cnt) = t2(i);
    It(cnt) = I0*exp(-(t2(i)-tmin)*R0/L0);
    cnt = cnt+1;
end
pp1 = interp1(t,It,'pchip','pp');
f1 = @(x) ppval(pp1,x);



fnum=fnum+1;
h(fnum)=figure;
plot(t.*1000.0,f1(t),'-b','linewidth',4)
 hold on
scatter(tmm,IR,10,'or','filled')
box on
xlabel('t (ms)','FontSize',18)
ylabel('I (A)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('Pure L/R','ANSYS Idump','Location','NorthEast')
set(gcf,'PaperPositionMode','auto')
% print(h(fnum),'-djpeg',['b1maglength', '-local'],'-r300')
print(h(fnum),'-djpeg','compare_LR','-r300')
hold off


%output Idump for this case
% tag = 'noLoss';
% dlmwrite(['IdumpLR_',tag,'.dat'],[tm,
    
%%% calculates MIITS for ansys data
pp4 = interp1(Iin(:,1),IR.*IR,'pchip','pp');
f4 = @(x) ppval(pp4,x);
MIITS_ansys = quad(f4,tmin,tmax)./(10^6);

Edump_ansys = MIITS_ansys*R0
Etot_L0 = L0_nofe*I0*I0/(2.0*10^6)


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Eddy Current Losses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hwedge = importdata('Hwedge.txt');
Hcollar = importdata('Hcollar.txt');
Hpole = importdata('Hpole.txt');

tmin = Hwedge(1,1);
tmax = Hwedge(length(Hwedge(:,1)),1);
t = tmin:1e-3:tmax

%%% calculates energy lost in wedges for ansys data
pp4 = interp1(Hwedge(:,1),Hwedge(:,2),'pchip','pp');
f4 = @(x) ppval(pp4,x);
E_wedge_ansys = quad(f4,tmin,tmax)


%%% calculates energy lost in collars for ansys data
pp5 = interp1(Hcollar(:,1),Hcollar(:,2),'pchip','pp');
f5 = @(x) ppval(pp5,x);
E_collar_ansys = quad(f5,tmin,tmax)

%%% calculates energy lost in poles for ansys data
pp6 = interp1(Hpole(:,1),Hpole(:,2),'pchip','pp');
f6 = @(x) ppval(pp6,x);
E_pole_ansys = quad(f6,tmin,tmax)


fnum=fnum+1;
h(fnum)=figure;
scatter(Hwedge(:,1)*1000,Hwedge(:,2)/1000.0,'ro','filled') 
hold on 
scatter(Hcollar(:,1)*1000,Hcollar(:,2)/1000.0,'bo','filled')  
scatter(Hpole(:,1)*1000,Hpole(:,2)/1000.0,'mo','filled') 
plot(t*1000,f4(t)/1000.0,'-r','linewidth',2)
plot(t*1000,f5(t)/1000.0,'-b','linewidth',2)
plot(t*1000,f6(t)/1000.0,'-m','linewidth',2)
box on
ylabel('Hgen (kW)','FontSize',18)
xlabel('time (ms)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('Wedges','Collars-nolam','Poles-nolam','Location','NorthEast')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','eddy_loss','-r300')
hold off


pdump = 100*Edump_ansys/Etot_L0 
pwedge = 100*E_wedge_ansys/(Etot_L0*10^6) 
pcollar = 100*E_collar_ansys/(Etot_L0*10^6) 
ppole = 100*E_pole_ansys/(Etot_L0*10^6) 

return




