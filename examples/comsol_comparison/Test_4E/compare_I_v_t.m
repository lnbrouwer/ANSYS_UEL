
clearvars 
close all
fnum = 0;

%convert circu res
run .\user102\write_trans_elem_res.m
datCirc = importdata('./user102/circu_res.txt');    %dlmwrite('circu_res.txt',[tm I1 V1 VR VS]);
I1 = abs(datCirc(:,2));
dat = importdata('./user102/Mag_v_t.txt');
tm1 = dat(:,1);



dat2 = importdata('../Test_4A/user102/circu_res.txt');
I2 = abs(dat2(:,2));
dat = importdata('../Test_4A/user102/Mag_v_t.txt');
tm2 = dat(:,1);

dat3 = importdata('../Test_4D/user102/circu_res.txt');
I3 = abs(dat3(:,2));
dat = importdata('../Test_4D/user102/Mag_v_t.txt');
tm3 = dat(:,1);



% fnum=fnum+1;
% h(fnum)=figure;
% scatter(tm2,I2,'filled')
% hold on 
% scatter(tm3,I3,'filled')
% scatter(tm1,I1,'filled')
% % scatter(sTHT0(:,1),sTHT0(:,2),10,'filled')
% box on
% hold on
% grid on
% % title(tag)
% xlabel('tm (s)','FontSize',18)
% ylabel('I (A)','FontSize',18)
% set(gca,'FontSize',16,'linewidth',2)
% set(h(fnum),'Position', [200 200 850 600])
% legend('ANSYS','ANSYS-IFCC','ANSYS-quenched','Location','NorthEast')
% set(gcf,'PaperPositionMode','auto')
% print(h(fnum),'-djpeg','compare_Ivt','-r300')
% hold off
% 
% 
% return 
% 


fnum=fnum+1;
h(fnum)=figure;
scatter(tm2,I2,'filled')
hold on 
scatter(tm3,I3,'filled')
scatter(tm1,I1,'filled')

% plot(tm2,I2,'--k','linewidth',1)
% hold on 
% plot(tm3,I3,'--k','linewidth',1)
% plot(tm1,I1,'--k','linewidth',1)


% scatter(sTHT0(:,1),sTHT0(:,2),10,'filled')
box on
hold on
% grid on
% title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('I (A)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('No Loss','IFCC only','IFCC w/Quench','Location','NorthEast')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','compare_Ivt','-r300')
hold off



%example of plot inside plot

t = linspace(0,2*pi,1000); % values in time
x = cos(t); % a sinusoid
perturbation = 0.1*exp((-(t-5*pi/4).^2)/.01).*sin(200*t); % a perturbation
signal = x+perturbation; % a signal to plot

% plot signal on new figure
figure,plot(t,x+perturbation)
xlabel('time'),ylabel('signal')

% create a new pair of axes inside current figure
axes('position',[.65 .175 .25 .25])
box on % put box around new pair of axes

indexOfInterest = (t < 11*pi/8) & (t > 9*pi/8); % range of t near perturbation

plot(t(indexOfInterest),signal(indexOfInterest)) % plot on new axes
axis tight



fnum=fnum+1;
h(fnum)=figure;
scatter(tm2,I2,'filled')
hold on 
scatter(tm3,I3,'filled')
scatter(tm1,I1,'filled')

% plot(tm2,I2,'--k','linewidth',1)
% hold on 
% plot(tm3,I3,'--k','linewidth',1)
% plot(tm1,I1,'--k','linewidth',1)

% create a new pair of axes inside current figure
axes('position',[.2 .4 .8 .8])
box on % put box around new pair of axes
indexOfInterest = (tm2 < 1e-3) & (tm2 >50e-3); % range of t near perturbation
% plot(tm2(indexOfInterest),I2(indexOfInterest)) % plot on new axes
plot(tm2(1:10),I2(1:10)) % plot on new axes
axis tight



% scatter(sTHT0(:,1),sTHT0(:,2),10,'filled')
box on
hold on
% grid on
% title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('I (A)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
legend('No Loss','IFCC only','IFCC w/Quench','Location','NorthEast')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','compare_Ivt_with_zoom','-r300')
hold off


return 


fnum=fnum+1;
h(fnum)=figure;
scatter(tm2,I2,'filled')
hold on 
scatter(tm3,I3,'filled')
scatter(tm1,I1,'filled')

plot(tm2,I2,'--k','linewidth',1)
hold on 
plot(tm3,I3,'--k','linewidth',1)
plot(tm1,I1,'--k','linewidth',1)


% scatter(sTHT0(:,1),sTHT0(:,2),10,'filled')
box on
hold on
% grid on
% title(tag)
xlabel('tm (s)','FontSize',18)
ylabel('I (A)','FontSize',18)
set(gca,'FontSize',16,'linewidth',2)
set(h(fnum),'Position', [200 200 850 600])
% legend('No Loss','IFCC only','IFCC w/Quench','Location','NorthEast')
set(gcf,'PaperPositionMode','auto')
print(h(fnum),'-djpeg','compare_noleg','-r300')
hold off







