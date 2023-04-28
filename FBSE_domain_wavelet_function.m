function ymw=FBSE_domain_wavelet_function(freq1, mb, ma, gamma, N)
w=freq1';
aw=w;
ymw=zeros(N,1);
wn=w(mb);
wm=w(ma);
an=1/(2*gamma*wn);
am=1/(2*gamma*wm);
pbn=(1+gamma)*wn;
mbn=(1-gamma)*wn;
pbm=(1+gamma)*wm;
mbm=(1-gamma)*wm;

for k=1:N
   if ((aw(k)>=pbn) && (aw(k)<=mbm))
       ymw(k)=1;
   elseif ((aw(k)>=mbm) && (aw(k)<=pbm))
       %ymw(k)=complex(cos(w(k)/2),-sin(w(k)/2))*cos(pi*EWT_beta(am*(aw(k)-mbm))/2);   
       ymw(k)=cos(pi*EWT_beta(am*(aw(k)-mbm))/2);   
   elseif ((aw(k)>=mbn) && (aw(k)<=pbn))
       %ymw(k)=complex(cos(w(k)/2),-sin(w(k)/2))*sin(pi*EWT_beta(an*(aw(k)-mbn))/2);
       ymw(k)=sin(pi*EWT_beta(an*(aw(k)-mbn))/2);
   end
end
end
