%est loss
rhocu = 1/(1e-8/pi);
area = 0.22500E-04;
I = 1e-6;

loss = rhocu*I*I/(area*area)
loss_sum = loss*area


%calc resisitivty for cern 
fcond = ns*pi*ds*ds/(4*aw*bw)
fsc = 0.4;
fcu = 1-fsc;
rhocu = (1e-8)/pi;
