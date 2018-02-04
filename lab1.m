clear
close all
%%
%Trial_1
mass=1;
trial;
[N,M]=size(datmat);
m=4;
for n=1:N
if datmat(n,m) ~= 0
firstsign=sign(datmat(n,m));
ns=n;
break;
end
end
peak=[];
neak=[];
numk=6;
for npk=1:numk
for n=ns+1:N
nsign=sign(datmat(n,m));
if nsign ~= firstsign
ns2=n;
break;
end
end
end
[ak,nk]=max(abs(datmat(ns:ns2,m)));
peak=[peak,ak*firstsign];
neak=[neak,ns+nk-1];
ns=ns2;
firstsign=nsign;
if ns>=N-1
end
h1=plot(datmat(:,1),datmat(:,m),'c-');
hold on
plot(neak-1,peak,'cx','linewidth',2);
%peak 1
x0=peak(1);
for npk=1:numk
for n=ns+1:N
nsign=sign(datmat(n,m));
if nsign ~= firstsign
ns3=n;
break;
end
end
end
[ak,nk]=max(abs(datmat(ns:ns3,m)));
peak=[peak,ak*firstsign];
neak=[neak,ns+nk-1];
ns=ns3;
firstsign=nsign;
if ns>=N-1
end
plot(neak-1,peak,'cx','linewidth',2);
ylabel('theta(t)')
xlabel('t')
title('System Response')
box on
grid on
axis on
%peak 2
x1=abs(peak(2));
for npk=1:numk
for n=ns+1:N
nsign=sign(datmat(n,m));
if nsign ~= firstsign
ns4=n;
break;
end
end
end
[ak,nk]=max(abs(datmat(ns:ns4,m)));
peak=[peak,ak*firstsign];
neak=[neak,ns+nk-1];
ns=ns4;
firstsign=nsign;
if ns>=N-1
end
plot(neak-1,peak,'cx','linewidth',2);
%peak 3
x2=peak(3);
T=neak(3)-neak(2);
delta=log(x1/x2);
Damp_ratio=real(1/(1+(2*pi/delta)^2)^0.5);
omega_d=2*pi/T/360;
omega_n=omega_d/(1-Damp_ratio^2)^0.5;
tau_d=delta/omega_n/Damp_ratio;
%% Trial 2
mass=1;
trial2;
[N,M]=size(datmat2);
m=4;
for n=1:N
if datmat2(n,m) ~= 0
firstsign=sign(datmat2(n,m));
ns=n;
break;
end
end
peak=[];
neak=[];
numk=6;
for npk=1:numk
for n=ns+1:N
nsign=sign(datmat2(n,m));
if nsign ~= firstsign
ns2=n;
break;
end
end
end
[ak,nk]=max(abs(datmat2(ns:ns2,m)));
peak=[peak,ak*firstsign];
neak=[neak,ns+nk-1];
ns=ns2;
firstsign=nsign;
if ns>=N-1
end
h2=plot(datmat(:,1),datmat2(:,m),'b-');
hold on
plot(neak-1,peak,'bx','linewidth',2);
%peak 1
x0=peak(1);
for npk=1:numk
for n=ns+1:N
nsign=sign(datmat2(n,m));
if nsign ~= firstsign
ns3=n;
break;
end
end
end
[ak,nk]=max(abs(datmat2(ns:ns3,m)));
peak=[peak,ak*firstsign];
neak=[neak,ns+nk-1];
ns=ns3;
firstsign=nsign;
if ns>=N-1
end
plot(neak-1,peak,'bx','linewidth',2);
ylabel('theta(t)')
xlabel('t')
title('System Response')
box on
grid on
axis on
%peak 2
x1=abs(peak(2));
for npk=1:numk
for n=ns+1:N
nsign=sign(datmat2(n,m));
if nsign ~= firstsign
ns4=n;
break;
end
end
end
[ak,nk]=max(abs(datmat2(ns:ns4,m)));
peak=[peak,ak*firstsign];
neak=[neak,ns+nk-1];
ns=ns4;
firstsign=nsign;
if ns>=N-1
end
plot(neak-1,peak,'bx','linewidth',2);
%peak 3
x2=peak(3);
T2=neak(3)-neak(2);
delta2=log(x1/x2);
Damp_ratio2=real(1/(1+(2*pi/delta2)^2)^0.5);
omega_d2=2*pi/T2/360;
omega_n2=omega_d2/(1-Damp_ratio2^2)^0.5;
tau_d2=delta2/omega_n2/Damp_ratio2;
%% Trial 3
mass=1;
trial3;
[N,M]=size(datmat3);
m=4;
for n=1:N
if datmat3(n,m) ~= 0
firstsign=sign(datmat3(n,m));
ns=n;
break;
end
end
peak=[];
neak=[];
numk=6;
for npk=1:numk
for n=ns+1:N
nsign=sign(datmat3(n,m));
if nsign ~= firstsign
ns2=n;
break;
end
end
end
[ak,nk]=max(abs(datmat3(ns:ns2,m)));
peak=[peak,ak*firstsign];
neak=[neak,ns+nk-1];
ns=ns2;
firstsign=nsign;
if ns>=N-1
end
h3=plot(datmat3(:,1),datmat3(:,m),'g-');
hold on
plot(neak-1,peak,'gx','linewidth',2);
%peak 1
x0=peak(1);
for npk=1:numk
for n=ns+1:N
nsign=sign(datmat3(n,m));
if nsign ~= firstsign
ns3=n;
break;
end
end
end
[ak,nk]=max(abs(datmat3(ns:ns3,m)));
peak=[peak,ak*firstsign];
neak=[neak,ns+nk-1];
ns=ns3;
firstsign=nsign;
if ns>=N-1
end
plot(neak-1,peak,'gx','linewidth',2);
ylabel('theta(t)')
xlabel('t')
title('System Response')
box on
grid on
axis on
%peak 2
x1=abs(peak(2));
for npk=1:numk
for n=ns+1:N
nsign=sign(datmat3(n,m));
if nsign ~= firstsign
ns4=n;
break;
end
end
end
[ak,nk]=max(abs(datmat3(ns:ns4,m)));
peak=[peak,ak*firstsign];
neak=[neak,ns+nk-1];
ns=ns4;
firstsign=nsign;
if ns>=N-1
end
plot(neak-1,peak,'gx','linewidth',2);
%peak 3
x2=peak(3);
T3=neak(3)-neak(2);
delta3=log(x1/x2);
Damp_ratio3=real(1/(1+(2*pi/delta3)^2)^0.5);
omega_d3=2*pi/T3/360;
omega_n3=omega_d3/(1-Damp_ratio3^2)^0.5;
tau_d3=delta3/omega_n3/Damp_ratio3;
%% Trial 4
mass=1;
trial4;
[N,M]=size(datmat4);
m=4;
for n=1:N
if datmat4(n,m) ~= 0
firstsign=sign(datmat4(n,m));
ns=n;
break;
end
end
peak=[];
neak=[];
numk=6;
for npk=1:numk
for n=ns+1:N
nsign=sign(datmat4(n,m));
if nsign ~= firstsign
ns2=n;
break;
end
end
end
[ak,nk]=max(abs(datmat4(ns:ns2,m)));
peak=[peak,ak*firstsign];
neak=[neak,ns+nk-1];
ns=ns2;
firstsign=nsign;
if ns>=N-1
end
h4=plot(datmat4(:,1),datmat4(:,m),'r-');
hold on
plot(neak-1,peak,'rx','linewidth',2);
%peak 1
x0=peak(1);
for npk=1:numk
for n=ns+1:N
nsign=sign(datmat4(n,m));
if nsign ~= firstsign
ns3=n;
break;
end
end
end
[ak,nk]=max(abs(datmat4(ns:ns3,m)));
peak=[peak,ak*firstsign];
neak=[neak,ns+nk-1];
ns=ns3;
firstsign=nsign;
if ns>=N-1
end
plot(neak-1,peak,'rx','linewidth',2);
ylabel('theta(t)')
xlabel('t')
title('System Response')
box on
grid on
axis on
%peak 2
x1=abs(peak(2));
for npk=1:numk
for n=ns+1:N
nsign=sign(datmat4(n,m));
if nsign ~= firstsign
ns4=n;
break;
end
end
end
[ak,nk]=max(abs(datmat4(ns:ns4,m)));
peak=[peak,ak*firstsign];
neak=[neak,ns+nk-1];
ns=ns4;
firstsign=nsign;
if ns>=N-1
end
plot(neak-1,peak,'rx','linewidth',2);
%peak 3
x2=peak(3);
T4=neak(3)-neak(2);
delta4=log(x1/x2);
Damp_ratio4=real(1/(1+(2*pi/delta4)^2)^0.5);
omega_d4=2*pi/T4/360;
omega_n4=omega_d4/(1-Damp_ratio4^2)^0.5;
tau_d4=delta4/omega_n4/Damp_ratio4;
%% Trial 5
mass=1;
trial5;
[N,M]=size(datmat5);
m=4;
for n=1:N
if datmat5(n,m) ~= 0
firstsign=sign(datmat5(n,m));
ns=n;
break;
end
end
peak=[];
neak=[];
numk=6;
for npk=1:numk
for n=ns+1:N
nsign=sign(datmat5(n,m));
if nsign ~= firstsign
ns2=n;
break;
end
end
end
[ak,nk]=max(abs(datmat5(ns:ns2,m)));
peak=[peak,ak*firstsign];
neak=[neak,ns+nk-1];
ns=ns2;
firstsign=nsign;
if ns>=N-1
end
h5=plot(datmat5(:,1),datmat5(:,m),'m-');
hold on
plot(neak-1,peak,'mx','linewidth',2);
%peak 1
x0=peak(1);
for npk=1:numk
for n=ns+1:N
nsign=sign(datmat5(n,m));
if nsign ~= firstsign
ns3=n;
break;
end
end
end
[ak,nk]=max(abs(datmat5(ns:ns3,m)));
peak=[peak,ak*firstsign];
neak=[neak,ns+nk-1];
ns=ns3;
firstsign=nsign;
if ns>=N-1
end
plot(neak-1,peak,'mx','linewidth',2);
ylabel('theta(t)')
xlabel('t')
title('System Response')
box on
grid on
axis on
%peak 2
x1=abs(peak(2));
for npk=1:numk
for n=ns+1:N
nsign=sign(datmat5(n,m));
if nsign ~= firstsign
ns4=n;
break;
end
end
end
[ak,nk]=max(abs(datmat5(ns:ns4,m)));
peak=[peak,ak*firstsign];
neak=[neak,ns+nk-1];
ns=ns4;
firstsign=nsign;
if ns>=N-1
end
plot(neak-1,peak,'mx','linewidth',2);
%peak 3
x2=peak(3);
T5=neak(3)-neak(2);
delta5=log(x1/x2);
Damp_ratio5=real(1/(1+(2*pi/delta5)^2)^0.5);
omega_d5=2*pi/T5/360;
omega_n5=omega_d5/(1-Damp_ratio5^2)^0.5;
tau_d5=delta5/omega_n5/Damp_ratio5;
%% Trial 6
mass=1;
trial6;
[N,M]=size(datmat6);
m=4;
for n=1:N
if datmat6(n,m) ~= 0
firstsign=sign(datmat6(n,m));
ns=n;
break;
end
end
peak=[];
neak=[];
numk=6;
for npk=1:numk
for n=ns+1:N
nsign=sign(datmat6(n,m));
if nsign ~= firstsign
ns2=n;
break;
end
end
end
[ak,nk]=max(abs(datmat6(ns:ns2,m)));
peak=[peak,ak*firstsign];
neak=[neak,ns+nk-1];
ns=ns2;
firstsign=nsign;
if ns>=N-1
end
h6=plot(datmat6(:,1),datmat6(:,m),'k-');
hold on
plot(neak-1,peak,'kx','linewidth',2);
%peak 1
x0=peak(1);
for npk=1:numk
for n=ns+1:N
nsign=sign(datmat6(n,m));
if nsign ~= firstsign
ns3=n;
break;
end
end
end
[ak,nk]=max(abs(datmat6(ns:ns3,m)));
peak=[peak,ak*firstsign];
neak=[neak,ns+nk-1];
ns=ns3;
firstsign=nsign;
if ns>=N-1
end
plot(neak-1,peak,'kx','linewidth',2);
ylabel('theta(t)')
xlabel('t')
title('System Response')
box on
grid on
axis on
%peak 2
x1=abs(peak(2));
for npk=1:numk
for n=ns+1:N
nsign=sign(datmat6(n,m));
if nsign ~= firstsign
ns4=n;
break;
end
end
end
[ak,nk]=max(abs(datmat6(ns:ns4,m)));
peak=[peak,ak*firstsign];
neak=[neak,ns+nk-1];
ns=ns4;
firstsign=nsign;
if ns>=N-1
end
plot(neak-1,peak,'kx','linewidth',2);
%peak 3
x2=peak(3);
T6=neak(3)-neak(2);
delta6=log(x1/x2);
Damp_ratio6=real(1/(1+(2*pi/delta6)^2)^0.5);
omega_d6=2*pi/T6/360;
omega_n6=omega_d6/(1-Damp_ratio6^2)^0.5;
tau_d6=delta6/omega_n6/Damp_ratio6;
%% Parameters
legend([h1 h2 h3 h4 h5 h6],'trial1','trial2','trial3','trial4','trial5','trial6','location','best')
Tmatrix=[T T2 T3 T4 T5 T6];
deltamatrix=[delta delta2 delta3 delta4 delta5 delta6];
Damp_ratiomatrix=[Damp_ratio Damp_ratio2 Damp_ratio3 Damp_ratio4 Damp_ratio5 Damp_ratio6];
omega_dmatrix=[omega_d omega_d2 omega_d3 omega_d4 omega_d5 omega_d6];
omega_nmatrix=[omega_n omega_n2 omega_n3 omega_n4 omega_n5 omega_n6];
taumatrix=[tau_d tau_d2 tau_d3 tau_d4 tau_d5 tau_d6];
fprintf('Tau sub d(trial 1,2,3,4,5,6 in seqence)'),disp(Tmatrix);
fprintf('Delta function(trial 1,2,3,4,5,6 in seqence)'),disp(deltamatrix);
fprintf('Damping ratio(trial 1,2,3,4,5,6 in seqence)'),disp(Damp_ratiomatrix);
fprintf('Omega sub d(trial 1,2,3,4,5,6 in seqence)'),disp(omega_dmatrix);
fprintf('Omega sub n(trial 1,2,3,4,5,6 in seqence)'),disp(omega_nmatrix);
ratio1=omega_n^2;
ratio2=omega_n2^2;
ratio3=omega_n3^2;
ratio4=omega_n4^2;
ratio5=omega_n5^2;
ratio6=omega_n6^2;
mc1=(0.5*ratio4-0.25*ratio1)/(ratio1-ratio4);
mc2=(0.5*ratio5-0.25*ratio2)/(ratio2-ratio5);
mc3=(0.5*ratio6-0.25*ratio3)/(ratio3-ratio6);
massmatrix=[mc1 mc2 mc3];
fprintf('Three masses obtained are'),disp(massmatrix);
massavg=(mc1+mc2+mc3)/3;
fprintf('the average mass is'),disp(massavg);
m1=massavg+.25;
m2=massavg+.25;
m3=massavg+.25;
m4=massavg+.5;
m5=massavg+.5;
m6=massavg+.5;
masstotalmatrix=[m1 m2 m3 m4 m5 m6];
fprintf('Six total masses obtained are'),disp(masstotalmatrix);
k1=(massavg+.25)*ratio1;
k2=(massavg+.25)*ratio2;
k3=(massavg+.25)*ratio3;
k4=(massavg+.5)*ratio4;
k5=(massavg+.5)*ratio5;
k6=(massavg+.5)*ratio6;
springmatrix=[k1 k2 k3 k4 k5 k6];
fprintf('six springs obtained are'),disp(springmatrix);
springavg=(k1+k2+k3+k4+k5+k6)/6;
fprintf('the average spring is'),disp(springavg);
c1=2*Damp_ratio*(springavg*(massavg+.25))^.5;
c2=2*Damp_ratio2*(springavg*(massavg+.25))^.5;
c3=2*Damp_ratio3*(springavg*(massavg+.25))^.5;
c4=2*Damp_ratio4*(springavg*(massavg+.5))^.5;
c5=2*Damp_ratio5*(springavg*(massavg+.5))^.5;
c6=2*Damp_ratio6*(springavg*(massavg+.5))^.5;
dampermatrix=[c1 c2 c3 c4 c5 c6];
fprintf('six damper obtained are'),disp(dampermatrix);
damperavg=(c1+c2+c3+c4+c5+c6)/6;
fprintf('the average damper is'),disp(damperavg);