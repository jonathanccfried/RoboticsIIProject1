clear

%% Initialize Variables
tfinal = 5;
step = 0.01;
q = [0.5;2.5;0]; 
qd = [0;2;0];
ez = [0;0;1];
rz = 0;
nr = 4; %3 to 5, number of range sensors
nb = 2; %1 to 3, number of bearing sensors
vcov = 1*[0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01]'; %Covariance of sensor noise
avcov = 1*[0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01]'; %Assumed Covariance of sensor noise
wcov = 1*[0.01 0.01]'; %Covariance of Input Noise
awcov = 1*[0.01 0.01]'; %Assumed Covariance of Input Noise
mapcov = 0.5*0.01; %Covariance of Landmark
amapcov = 0.5*0.01;
k1 = 5;
k2 = 3; %Controller Variables
P = 0.15*eye(3+2*(nb+nr));
Sigma = 0.15*eye(3+2*(nb+nr));
Sigma(4:end,4:end) = amapcov*Sigma(4:end,4:end);
% Range Sensors, from 3 to 5

pL = [0 0 0;
      10 10 0;
      0 10 0;
      10 0 0;
      5 5 0]';

% Bearing Sensors, from 1 to 3

pB = [10 10 0
      10 5   0;
      5 10  0]';

LStateTrue = [];

for i=1:nr
    LStateTrue = [LStateTrue;pL(1:2,i)];
end

for i=1:nb
    LStateTrue = [LStateTrue;pB(1:2,i)];
end
%% Initialize lists

tlist = [];
qdlist = [];
dqdlist = [];
qlist = [];
LStatelist=[];
qhatKFlist = [];
ulist = [];
ylist = [];
utruelist = [];
Plist = [];
Sigmalist = [];

%% Initial Computations

LState = LStateTrue + randn(size(LStateTrue))*mapcov;

%% Main Loop
for t=0:step:tfinal
%% Store time

tlist = [tlist t];
    
    %% Define Trajectory
 if t < 2.5

    vd = 3;
    wd = pi/5;

    dqd = [vd*cos(qd(3));
           vd*sin(qd(3));
         wd];
%    dqd = [vd;wd];
 else
     vd = 2;
     wd = -pi/2;
    dqd = [vd*cos(qd(3));
           vd*sin(qd(3));
         wd];

     % 
% 
%    qd = [5 + vd*cos(wd*t);
%         5 + vd*sin(wd*t)
%         pi/2 + wd*t];
%     dqd = [vd;wd];
% 
 end

%% Calculate Sensor Data
pR = [q(1:2);0];
theta = q(3);
d = [0 0 0 0 0 0 0 0]';

%range sensors

for i=1:size(pL,2)
    d(i)=norm(pR-pL(:,i));
end

%bearing sensors

for i=1:size(pB,2)
    d(5+i)=atan2(pB(2,i)-q(2),pB(1,i)-q(1))-theta;
end
%% Disturb Sensor Data
y = d + randn(size(d)).*vcov;
%% Calculate Estimation - Mean Approach - at t0 to initialize the kalman filter

if t==0

% use two range sensors intersecting with the x-y plane
qhatexp = [0;0;0];
count =0;
pLhat = zeros(3,nr);
pBhat = zeros(3,nb);
for i=1:nr
pLhat(1:2,i) = LState(2*i-1:2*i);
end
for i=1:nb
pBhat(1:2,i) = LState(2*nr+2*i-1:2*nr+2*i);
end
for i=1:nr-2
for j=i+1:nr-1
for k=j+1:nr
kvec=(pL(:,i)-pL(:,j))/norm(pL(:,i)-pL(:,j));
d1=sqrt(y(i)^2-(pL(:,i)'*ez-rz)^2);
d2=sqrt(y(j)^2-(pL(:,j)'*ez-rz)^2);
phi=subprob3(ez,d1*kvec,pL(:,j)-pL(:,i),d2);
pp1=pL(:,i)-pL(:,i)'*ez*ez;
qest=[pp1+rot(ez,phi(1))*d1*kvec+rz*ez pp1+rot(ez,-phi(1))*d1*kvec+rz*ez];

opt = abs([norm(qest(:,1) - pL(:,k))-y(k), norm(qest(:,2)- pL(:,k)-y(k))]);
[dum,ind] = min(opt);


% use three balls intersection to disambiguate
% qest1=threeballs(pL(:,i),pL(:,j),pL(:,k),y(i),y(j),y(k));
% 
% [dum,ind]=min(vecnorm(qest(1:2,:)-mean(qest1(1:2,:),2)));
qhatexp= qhatexp + qest(:,ind);
count = count+1;
end
end
end

qhat=qhatexp(1:2)/count;
if isnan(qhat)   %If the balls don't intesect due to noise, just repeat last position
    qhat = qhatlist(1:2,end-1);
end

% calculate angle with bearing sensor

qhatexp = 0;
for i=1:nb
qhatexp=qhatexp + atan2(pB(2,i)-qhat(2,1),pB(1,i)-qhat(1,1))-y(5+i);
end
qhat(3,1) = qhatexp/nb;
end

%% Calculate Estimation - Kalman Filter Initial

if t==0
    qhatKF = qhat;
    KalmanState = [qhatKF;LState];
end

%% Calculate Error
e = qd - qhatKF;
%e = qd - q;
e_xy = e(1:2);
etheta = e(3);
%% Calculate Control


vl = vd;
wl = wd;
Rn = [cos(q(3)) sin(q(3))]';
Rl = [-sin(q(3)) cos(q(3))]';

v = k1*e_xy'*Rn + vl*cos(etheta);
w = vl*e_xy'*Rl + k2*sin(etheta) + wl;

u = [v;w];

%% Calculate Estimation - Kalman Filter Final

    qkalman = qhatKF + [cos(qhatKF(3)) 0; sin(qhatKF(3)) 0;0 1]*u*step;
    A=[1 0 -u(1)*sin(qhatKF(3)) zeros(1,2*(nb+nr));0 1 u(1)*cos(qhatKF(3)) zeros(1,2*(nb+nr)); 0 0 1 zeros(1,2*(nb+nr));zeros(2*(nb+nr),3) eye(size(LState*LState'))];
    B=[cos(qhatKF(3)) 0 ; sin(qhatKF(3)) 0 ; 0 1 ;zeros(2*(nb+nr),2)];
    M = diag([avcov(1:nr);avcov(6:5+nb)]);
    N = diag(awcov);
    yhat= zeros(nr+nb,1);
    gradienth = zeros(nr+nb,3+2*(nb+nr));
    for j=1:nr
    yhat(j) = norm(qkalman(1:2) - LState(2*j-1:2*j));
    gradienth(j,1:3) = [-(LState(2*j-1:2*j) - qkalman(1:2))'/norm((LState(2*j-1:2*j) - qkalman(1:2))),0];
    gradienth(j,3+2*j-1:3+2*j) = [(LState(2*j-1:2*j) - qkalman(1:2))'/norm(LState(2*j-1:2*j) - qkalman(1:2))];
    end
    for j=1:nb
    yhat(nr+j) = atan2(LState(2*j+2*nr)-qkalman(2),LState(2*j+2*nr-1)-qkalman(1))-qkalman(3);
    gradienth(nr+j,1:3)= [(LState(2*j+2*nr)-qkalman(2))/(norm(LState(2*j+2*nr-1:2*j+2*nr)-qkalman(1:2))^2),-(LState(2*j+2*nr-1)-qkalman(1))/(norm(LState(2*j+2*nr-1:2*j+2*nr)-qkalman(1:2))^2),-1];
    gradienth(nr+j,3+2*nr+2*j-1:3+2*nr+2*j)= [-(LState(2*j+2*nr)-qkalman(2))/(norm(LState(2*j+2*nr-1:2*j+2*nr)-qkalman(1:2))^2),(LState(2*j+2*nr-1)-qkalman(1))/(norm(LState(2*j+2*nr-1:2*j+2*nr)-qkalman(1:2))^2)];
    end
    C=gradienth;
    K = P*C'*inv(M);
    KalmanState = [qkalman;LState] + K*([y(1:nr);y(6:5+nb)] - yhat);
    qhatKF = KalmanState(1:3);
    LState = KalmanState(4:end);
    P=inv(inv(Sigma)+C'*inv(M)*C);
    Sigma=A*P*A'+B*N*B';
%% Disturb Control
utrue = u + randn(size(u)).*wcov;
%% Store Values
qdlist = [qdlist qd];
dqdlist = [dqdlist dqd];
qlist = [qlist q];
ylist = [ylist y];
ulist = [ulist u];
utruelist = [utruelist utrue];
LStatelist=[LStatelist LState];
qhatKFlist = [qhatKFlist qhatKF];
Plist = [Plist P];
Sigmalist = [Sigmalist Sigma];
%% Update States
q = q + [cos(q(3)) 0; sin(q(3)) 0;0 1]*utrue*step;
qd = qd + dqd*step;
%% Close

end