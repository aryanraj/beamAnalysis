%{
STRUCTURE of data.dat

0(length)         <length>
1(conc load)      <distance> <force>
2(conc moment)    <distance> <moment>
3(UDL)            <distance 1> <distance 2> <force>
4(LDL)            <distance 1> <distance 2> <force 1> <force 2>
5(hindge support) <distance>
6(fixed support)  <distance>
7(moment release) <distance>
8(shear release)  <distance>
9(EI)             <EI value>

%}

%t = cputime;

nParts = 1000;
M = zeros(nParts,1);
F = zeros(nParts,1);
EI = ones(nParts,1);
element = zeros(nParts,1);

%Take INPUT
data = importdata('data.dat');
for i=1:size(data)(1)
  switch data(i)
    case 0 %length
      L = data(i,2);
      L_ = L/(nParts-1);
    case 1 %concentrated load
      d = floor(data(i,2)/L_+1);
      f = data(i,3);
      F(d,1) = F(d,1)+f;
    case 2 %concntrated moment
      d = floor(data(i,2)/L_+1);
      m = data(i,3);
      M(d,1) = M(d,1)+m;
    case 3 %UDL
      d1 = floor(data(i,2)/L_+1);
      d2 = floor(data(i,3)/L_+1);
      f = data(i,4);
      F(d1:d2,1) = F(d1:d2,1) + L_*f;
    case 4 %LDL
      d1 = floor(data(i,2)/L_+1);
      d2 = floor(data(i,3)/L_+1);
      f1 = data(i,4);
      f2 = data(i,5);
      F(d1:d2,1) = F(d1:d2,1) + (f1:(f2-f1)/(d2-d1):f2)'*L_;
    case 5 %hindge support
      d = floor(data(i,2)/L_+1);
      element(d) = 1;
    case 6 %fixed support
      d = floor(data(i,2)/L_+1);
      element(d) = 2;
    case 7 %moment release
      d = floor(data(i,2)/L_+1);
      element(d) = 3;
    case 8 %shear release
      d = floor(data(i,2)/L_+1);
      element(d) = 4;
    case 9 %EI
      f = data(i,2);
      EI *= f;
  end
end

%Pre-processing of Data
SF = zeros(nParts,1);
BM = zeros(nParts,1);
Angle = zeros(nParts,1);
v = zeros(nParts,1);

%things in start
t = zeros(nParts,1);
t(1) = 1;
switch element(1)
  case 0
    Angle = [Angle t];
    v = [v t];
  case 1
    Angle = [Angle t];
    SF = [SF t];
  case 2
    SF = [SF t];
    BM = [BM t];
end

%things in between
for i=2:(nParts-1)
  switch element(i)
    case 1
      t = zeros(nParts,1);
      t(i) = 1;
      SF = [SF t];
    case 2
      t = zeros(nParts,1);
      t(i) = 1;
      SF = [SF t];
      SM = [BM t];
    case 3
      t = zeros(nParts,1);
      t(i) = 1;
      Angle = [Angle t];
    case 4
      t = zeros(nParts,1);
      t(i) = 1;
      v = [v t];
  end
end

SF(1) += F(1);
BM(1) += M(1);

BM = [BM(:,1) zeros(nParts,size(SF)(2)-1) BM(:,2:end)];
Angle = [Angle(:,1) zeros(nParts,size(BM)(2)-1) Angle(:,2:end)];
v = [v(:,1) zeros(nParts,size(Angle)(2)-1) v(:,2:end)];

%Form Equations
size_SF = size(SF)(2);
size_BM = size(BM)(2);
size_Angle = size(Angle)(2);
size_v = size(v)(2);

for i=2:nParts
  SF(i,:) += SF(i-1,:);
  SF(i) -= F(i);
  BM(i,:) += BM(i-1,:);
  BM(i,1:size_SF) += SF(i,:)*L_;
  BM(i) += M(i);
  Angle(i,:) += Angle(i-1,:);
  Angle(i,1:size_BM) += BM(i,:)/EI(i)*L_;
  v(i,:) += v(i-1,:);
  v(i,1:size_Angle) += Angle(i,:)*L_;
end

static_eqn_A = [];
static_eqn_B = [];
kinematic_eqn_A = [];
kinematic_eqn_B = [];
overall_eqn_A = [];
overall_eqn_B = [];

%in between
for i=2:(nParts-1)
  switch element(i)
    case 1
      %v=0
      t = zeros(1,size_v-1);
      t(1,1:(size_v-1)) += v(i,2:end);
      kinematic_eqn_A = [kinematic_eqn_A;t];
      kinematic_eqn_B = [kinematic_eqn_B;(0-v(i))];
    case 2
      %v=0
      t = zeros(1,size_v-1);
      t(1,1:(size_v-1)) += v(i,2:end);
      kinematic_eqn_A = [kinematic_eqn_A;t];
      kinematic_eqn_B = [kinematic_eqn_B;(0-v(i))]; 
      %Angle=0
      t = zeros(1,size_v-1);
      t(1,1:(size_Angle-1)) += Angle(i,2:end);
      kinematic_eqn_A = [kinematic_eqn_A;t];
      kinematic_eqn_B = [kinematic_eqn_B;(0-Angle(i))];
    case 3
      %BM=0
      t = zeros(1,size_BM-1);
      t(1,1:(size_BM-1)) += BM(i,2:end);
      static_eqn_A = [static_eqn_A;t];
      static_eqn_B = [static_eqn_B;(0-BM(i))];
    case 4
      %SF=0
      t = zeros(1,size_BM-1);
      t(1,1:(size_SF-1)) += SF(i,2:end);
      static_eqn_A = [static_eqn_A;t];
      static_eqn_B = [static_eqn_B;(0-SF(i))];
  end
end

%at End
switch element(nParts)
  case 0
    %SF=0
    t = zeros(1,size_BM-1);
    t(1,1:(size_SF-1)) += SF(i,2:end);
    static_eqn_A = [static_eqn_A;t];
    static_eqn_B = [static_eqn_B;(0-SF(i))];
    %BM=0
    t = zeros(1,size_BM-1);
    t(1,1:(size_BM-1)) += BM(i,2:end);
    static_eqn_A = [static_eqn_A;t];
    static_eqn_B = [static_eqn_B;(0-BM(i))];
  case 1
    %v=0
    t = zeros(1,size_v-1);
    t(1,1:(size_v-1)) += v(i,2:end);
    kinematic_eqn_A = [kinematic_eqn_A;t];
    kinematic_eqn_B = [kinematic_eqn_B;(0-v(i))];
    %BM=0
    t = zeros(1,size_BM-1);
    t(1,1:(size_BM-1)) += BM(i,2:end);
    static_eqn_A = [static_eqn_A;t];
    static_eqn_B = [static_eqn_B;(0-BM(i))];
  case 2
    %v=0
    t = zeros(1,size_v-1);
    t(1,1:(size_v-1)) += v(i,2:end);
    kinematic_eqn_A = [kinematic_eqn_A;t];
    kinematic_eqn_B = [kinematic_eqn_B;(0-v(i))]; 
    %Angle=0
    t = zeros(1,size_v-1);
    t(1,1:(size_Angle-1)) += Angle(i,2:end);
    kinematic_eqn_A = [kinematic_eqn_A;t];
    kinematic_eqn_B = [kinematic_eqn_B;(0-Angle(i))];
end

% static_eqn_A
% static_eqn_B
% kinematic_eqn_A
% kinematic_eqn_B

overall_eqn_A = zeros(size_v-1);
overall_eqn_A(1:size(static_eqn_A)(1),1:size(static_eqn_A)(2)) = static_eqn_A;
overall_eqn_A((size(static_eqn_A)(1)+1):end,:) = kinematic_eqn_A;
overall_eqn_B = [static_eqn_B;kinematic_eqn_B];

%overall_eqn_A
%overall_eqn_B

%calculating
det_overall = det(overall_eqn_A);

if issquare(static_eqn_A) == 0
  if det_overall == 0
    disp('Unstable Structure');
  else
    disp('Indeterminate Structure');
  end
end

overall_sol=[];

if det_overall ~= 0
  overall_sol = (inv(overall_eqn_A) * overall_eqn_B)';
end

final_SF = zeros(nParts,1);
final_BM = zeros(nParts,1);
final_Angle = zeros(nParts,1);
final_v = zeros(nParts,1);

%post calcuation
for i=1:nParts
  final_SF(i) = sum([SF(i) (SF(i,2:end).*overall_sol(1:(size_SF-1)))],2);
  final_BM(i) = sum([BM(i) (BM(i,2:end).*overall_sol(1:(size_BM-1)))],2);
  final_Angle(i) = sum([Angle(i) (Angle(i,2:end).*overall_sol(1:(size_Angle-1)))],2);
  final_v(i) = sum([v(i) (v(i,2:end).*overall_sol(1:(size_v-1)))],2);
end

%Plotting
figure
plot((1:nParts)*L_,final_SF), title('Shear Force')
figure
plot((1:nParts)*L_,final_BM), title('Bending Moment')
figure
plot((1:nParts)*L_,final_Angle), title('Angle')
figure
plot((1:nParts)*L_,final_v), title('Deflection')

%e = cputime-t

