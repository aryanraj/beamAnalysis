nParts = 10000;
M = zeros(nParts,1);
F = zeros(nParts,1);
cond_dist = [];
cond_moment = [];

%Take data
data = importdata('data.dat');

for i=1:size(data)(1)
  switch data(i,1)
    case 0
      L = data(i,2);
      L_ = L/nParts;
    case 1
      d = data(i,2);
      f = data(i,3);
      F(int16(d/L_),1) = F(int16(d/L_),1)+f;
    case 2
      d = data(i,2);
      m = data(i,3);
      M(int16(d/L_),1) = M(int16(d/L_),1)+m;
    case 3
      d1 = data(i,2);
      d2 = data(i,3);
      f = data(i,4);
      F(int16(d1/L_):int16(d2/L_),1) = F(int16(d1/L_):int16(d2/L_),1) + L_*f;
    case 4
      d1 = data(i,2);
      d2 = data(i,3);
      f1 = data(i,4);
      f2 = data(i,5);
      F(int16(d1/L_):int16(d2/L_),1) = F(int16(d1/L_):int16(d2/L_),1) + (f1:(f2-f1)/((f2-f1)/L_):f2)*L_;
    case 5
      d = data(i,2);
      f = zeros(nParts);
      f(int16(d/L_)) = 1;
      F = [F f];
    case 6
      
    case 7
      
    case 8
      
  end
end

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

n=input('Total number of reactions : ');
for i=1:n
  d=input('dist : ');
  F=[F zeros(nParts,1)];
  F(i+1,int16(d/L_))=-1;
end

BM = zeros(nParts,BM_vars+1);
EI = ones(nParts,SF_vars+1);
Angle = zeros(nParts,Angle_vars+1);
v = zeros(nParts,v_vars+1);

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
