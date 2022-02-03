clear

%% Initialize Variables
tfinal = 5;
step = 0.01;
q = [5;7.5;0]; 
qd = [5;8;0];
ez = [0;0;1];
rz = 0;
nr = 3; %3 to 5, number of range sensors
nb = 1; %1 to 3, number of bearing sensors
EM = 3; %1 - Mean, 2 - NL Least Squares, 3 -  Kalman Filter
vcov = 5*[0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01]'; %Covariance of sensor noise
avcov = 5*[0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01]'; %Assumed Covariance of sensor noise
wcov = 5*[0.01 0.01]'; %Covariance of Input Noise
awcov = 5*[0.01 0.01]'; %Assumed Covariance of Input Noise
mapcov = 0.01; %Covariance of Landmark
amapcov = 0.01;
NLN = 40;
k1 = 5;
k2 = 3; %Controller Variables
P = eye(2*(nb+nr));
Sigma = eye(2*(nb+nr));

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

    vd = 2;
    wd = pi/10;

    dqd = [vd*cos(qd(3));
           vd*sin(qd(3));
         wd];
%    dqd = [vd;wd];
 else
     vd = .5;
     wd = -pi/6;
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



%% Calculate Error
e = qd - q;
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
%% Disturb Control
utrue = u + randn(size(u)).*wcov;

%% Calculate Mapping Estimation - Kalman Filter Final

    A=eye(size(LState*LState'));
    M = diag([avcov(1:nr);avcov(6:5+nb)]);
 %   for i=1:NLN
    yhat= zeros(nr+nb,1);
    gradienth = zeros((nr+nb),2*(nr+nb));
    for j=1:nr
    yhat(j) = norm(q(1:2) - LState(2*j-1:2*j));
    gradienth(j,2*j-1:2*j) = [(LState(2*j-1:2*j) - q(1:2))'/norm(LState(2*j-1:2*j) - q(1:2))];
    end
    for j=1:nb
        yhat(nr+j) = atan2(LState(2*j+2*nr)-q(2),LState(2*j+2*nr-1)-q(1))-q(3);
        gradienth(nr+j,2*nr+2*j-1:2*nr+2*j)= [-(LState(2*j+2*nr)-q(2))/(norm(LState(2*j+2*nr-1:2*j+2*nr)-q(1:2))^2),(LState(2*j+2*nr-1)-q(1))/(norm(LState(2*j+2*nr-1:2*j+2*nr)-q(1:2))^2)];
    end
  %  end
    C=gradienth;
    K = P*C'*inv(M);
    LState = LState + K*([y(1:nr);y(6:5+nb)] - yhat);
    P=inv(inv(Sigma)+C'*inv(M)*C);
    Sigma=A*P*A';
%% Store Values
qdlist = [qdlist qd];
dqdlist = [dqdlist dqd];
qlist = [qlist q];
ylist = [ylist y];
utruelist = [utruelist utrue];
LStatelist = [LStatelist LState];
Plist = [Plist P];
Sigmalist = [Sigmalist Sigma];
%% Update States
q = q + [cos(q(3)) 0; sin(q(3)) 0;0 1]*utrue*step;
qd = q + dqd*step;
%% Close

end