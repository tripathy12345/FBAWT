function yms=FBSE_domain_scaling_function(freq1, boundaries, gamma, N)
w=freq1';
w1=w(boundaries(1));
aw=w;
yms=zeros(N,1);

an=1/(2*gamma*w1);
pbn=(1+gamma)*w1;
mbn=(1-gamma)*w1;
for k=1:N
   if (aw(k)<=mbn)
       yms(k)=1;
   elseif ((aw(k)>=mbn) && (aw(k)<=pbn))
       yms(k)=cos(pi*EWT_beta(an*(aw(k)-mbn))/2);
   end
end
end
