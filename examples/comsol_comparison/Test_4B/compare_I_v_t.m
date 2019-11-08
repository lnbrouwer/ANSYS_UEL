

clearvars 
close all
fnum = 0;

%convert circu res
run .\user102\write_trans_elem_res.m
datCirc = importdata('./user102/circu_res.txt');    %dlmwrite('circu_res.txt',[tm I1 V1 VR VS]);
I1 = abs(datCirc(:,2));
dat = importdata('./user102/Mag_v_t.txt');
tm1 = tm;



dat2 = importdata('../Test_4A/user102/circu_res.txt');
I2 = abs(dat2(:,2));
dat = importdata('../Test_4A/user102/Mag_v_t.txt');
tm2 = tm;


fnum=fnum+1;
h(fnum)=figure;
scatter(tm1,I1,'filled')
hold on 
scatter(tm2,I2,'filled')
% scatter(sTHT0(:,1),sTHT0(:,2),10,'filled')
box on
hold on
grid on
% title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('I (A)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('ANSYS','ANSYS-IFCC','Location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','compare_Ivt','-r300')
hold off













