clearvars 
close all
fnum = 0;


tag = 'Test3B';


dat2 = importdata('../Test_3Ba/user102/Therm_v_t.txt');
T1 = dat2(:,2);
C1 = dat2(:,5);

dat2 = importdata('../Test_3Bb/user102/Therm_v_t.txt');
T2 = dat2(:,2);
C2 = dat2(:,5);


dat2 = importdata('./user102/Therm_v_t.txt');
T3 = dat2(:,2);
C3 = dat2(:,5);



dat2 = importdata('..\..\CERN_res\3Ba\1e6\CvSc_avg_HT0_5Ba.csv');
dat22 = importdata('..\..\CERN_res\3Ba\1e6\T_avg_HT0_5Ba.csv');

T4 = dat22(:,2);
C4 = dat2(:,2);

dat2 = importdata('..\..\CERN_res\3Bb\1e6\CvSc_avg_HT0_3Bb.csv');
dat22 = importdata('..\..\CERN_res\3Bb\1e6\T_avg_HT0_3Bb.csv');
T5 = dat22(:,2);
C5 = dat2(:,2);

dat2 = importdata('..\..\CERN_res\3Bc\1e6\CvSc_avg_HT0_3Bc.csv');
dat22 = importdata('..\..\CERN_res\3Bc\1e6\T_avg_HT0_3Bc.csv');
T6 = dat22(:,2);
C6 = dat2(:,2);




fnum=fnum+1;
h(fnum)=figure;
scatter(T1,C1,'filled')
hold on 
scatter(T2,C2,'filled')
scatter(T3,C3,'filled')
plot(T4,C4)
plot(T5,C5)
plot(T6,C6)

box on
hold on
grid on
title([tag,' compare 0,5,10 T'])
axis([12 20 0 2e5])
xlabel('temp (K)','FontSize',18)
ylabel('ave cvNb3Sn HT0 (J/kg*m^3)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('B=0 ANSYS','B=5 ANSYS','B=10 ANSYS','B=0 COMSOL','B=5 COMSOL','B=10 COMSOL','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','check_B_var','-r300')
hold off







