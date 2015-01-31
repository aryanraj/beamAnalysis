nParts = 10000;
M = zeros(1,nParts);
F = zeros(1,nParts);
SF = zeros(1,nParts);
BM = zeros(1,nParts);
EI = ones(1,nParts);
Angle = zeros(1,nParts);
v = zeros(1,nParts);

%Take data
L = input('Length : ');
L_ = L/nParts;

EI = input('EI : ')*ones(1,nParts);


n = input('Total conc. loads : ');
for i=1:n
  d = input('dist : ');
  F(int16(d/L_)) = F(int16(d/L_))+input('load : ');
end

n = input('Total conc. moments : ');
for i=1:n
  d = input('dist : ');
  M(int16(d/L_)) = M(int16(d/L_))+input('moment : ');
end

n = input('Total uniform distributed load : ');
for i=1:n
  d1=input('start : ');
  d2=input('end : ');
  F(int16(d1/L_):int16(d2/L_)) = F(int16(d1/L_):int16(d2/L_)) + L_*input('load per length:');
end

n = input('Total linear distributed load : ');
for i=1:n
  d1=input('start : ');
  d2=input('end : ');
  l1=input('start load : ');
  l2=input('end load : ');
  %disp((l1:(l2-l1)/((d2-d1)/L_):l2));
  F(int16(d1/L_):int16(d2/L_)) = F(int16(d1/L_):int16(d2/L_)) + (l1:(l2-l1)/((d2-d1)/L_):l2)*L_;
end

t = cputime;

s=0;
for i=1:nParts
  s=s+F(i)*i*L_+M(i);
end
F(nParts)=F(nParts)-s/L;
F(1)=F(1)-sum(F);

SF(1)=-F(1);
BM(1)=M(1);
for i=2:nParts
  SF(i)=SF(i-1)-F(i);
  BM(i)=BM(i-1)+SF(i)*L_+M(i);
  Angle(i)=Angle(i-1)+BM(i)/EI(i)*L_;
  v(i)=v(i-1)+Angle(i)*L_;
end
SF(1)=0;

%correction
corr = v(nParts)/L;
for i=1:nParts
  Angle(i)=Angle(i)-corr;
  v(i)=v(i)-corr*i*L_;
end

e = cputime-t
