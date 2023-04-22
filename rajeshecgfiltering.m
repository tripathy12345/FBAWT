function zz=rajeshecgfiltering(xx,fs)
%%%%%%%%%butterworth HPF design%%%%%%
[b,a]=butter(1,(0.5)/(fs/2), 'high');
yy=filtfilt(b,a,xx);
%%%%%%%%%butterworth HPF design%%%%%%
[b1,a1]=butter(1,(45)/(fs/2), 'low');
zz=filtfilt(b1,a1,yy);
%  yy=yy/max(abs(yy));
end