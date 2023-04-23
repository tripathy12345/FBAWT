clc;
clear all;
close all;
load ecg_signal.mat;  
Fs=500;
zf=ecg_signal;
zf=ecgfiltering(zf,Fs);
f=zf;
figure,
plot(zf)
figure
N=length(f);
nb=(1:N);
f=f';
MM=N;
%%%%%FBSE coefficents evaluation in an iterative way%%%%%%
if exist('alfa') == 0
    x=2;
    alfa=zeros(1,MM);
    for i=1:MM
        ex=1;
        while abs(ex)>.00001
            ex=-besselj(0,x)/besselj(1,x);
            x=x-ex;
        end
        alfa(i)=x;
        fprintf('Root # %g  = %8.5f ex = %9.6f \n',i,x,ex)
        x=x+pi;
    end
end
a=N;
%%%%FBSE Coefficents%%%%%%%%
for m1=1:MM
    a3(m1)=(2/(a^2*(besselj(1,alfa(m1))).^2))*sum(nb.*f.*besselj(0,alfa(m1)/a*nb));
end
%%%%root to frequency conversion%%%%%%%%
freq1=(alfa*Fs)/(2*pi*length(f));
order=2;
framelen=(N/(0.1*Fs))+1;
%%%%Instantaneous energy in FBSE domain%%%%%
eng=((besselj(1,alfa(m1)).^2).*(length(f)^2).*(a3.^2))/2;
a3new=sgolayfilt(eng,order,framelen);
% plot(freq1, a3new)
Nb=100;
bound = LocalMax(a3new',Nb);
fbound=freq1(bound);
boundaries=bound;
Npic=length(boundaries);
% We compute gamma accordingly to the theory
gamma=1;
for k=1:Npic-1
    r=(boundaries(k+1)-boundaries(k))/(boundaries(k+1)+boundaries(k));
    if r<gamma 
       gamma=r;
    end
end
mfb=cell(Npic+1,1);
%%%%FBSE domain scaling function (No concept of FFT here like EWT)%%%
mfb{1}=FBSE_domain_scaling_function(freq1, boundaries, gamma, N);
%%%%FBSE domain wavelet function (No concept of FFT here like EWT)%%%
for k=1:Npic-1
   mfb{k+1}=FBSE_domain_wavelet_function(freq1,boundaries(k),boundaries(k+1),gamma,N); 
end
mfb{Npic+1}=FBSE_domain_wavelet_function(freq1,boundaries(Npic),length(a3),gamma,N); 
%%%%%show_boundaries%%%%%%%%%5
plot(freq1,a3);
hold on
for i=1:Npic
     line([freq1(boundaries(i)) freq1(boundaries(i))],[0 max(abs(a3))],'LineWidth',2,'LineStyle','--','Color',[1 0 0]);
end
figure,
for i=1:size(mfb)-1
plot(freq1,mfb{i,1}, 'k', 'LineWidth', 1)
hold on
end
xlim([0 45])
ylim([0 2])
title('FB domain-adaptive wavelet filter bank')
%%%%FBSE domain mode reconstruction (No concept of IFFT here as in EWT)%%%
for m1=1:MM
D(:,m1)=besselj(0,alfa(m1)/a*nb);
end
fbseewt=cell(length(mfb),1);
for k=1:length(mfb)
    fbseewt{k}=((mfb{k})'.*a3);
    modes{k}=fbseewt{k}*D;
end
figure
subplot(811)
plot(modes{1});
subplot(812)
plot(modes{2});
subplot(813)
plot(modes{3});
subplot(814)
plot(modes{4});
subplot(815)
plot(modes{5});
subplot(816)
plot(modes{6});
subplot(817)
plot(modes{7});
subplot(818)
plot(modes{8});
