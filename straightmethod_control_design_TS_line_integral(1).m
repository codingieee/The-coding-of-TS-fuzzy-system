clear
clc
% 
yalmip('clear')
mu = 0.5;
h = 1373;

check_lmi(mu,h)

function feasibli = check_lmi(mu,h)

%% system dimension

rule = 2;

A = cell(1, rule);
Ad = cell(1, rule);

% %%%%%%%%%%%%%%%%%%%%%% 例子1 
% % Parameters
% L     = 5.5;
% l     = 2.8;
% vbar  = -1.0;
% tbar  = 2.0;
% t0    = 0.5;
% rho   = 0.7;
% kappa = 1e-2;
% 
% % System matrices
% A{1} = [ -rho*vbar*tbar/(L*t0),      0,        0;
%         rho*vbar*tbar/(L*t0),      0,        0;
%         rho*vbar^2*tbar^2/(2*L*t0), vbar*tbar/t0, 0 ];
% 
% Ad{1} = [ -(1-rho)*vbar*tbar/(L*t0),      0,        0;
%          (1-rho)*vbar*tbar/(L*t0),      0,        0;
%          (1-rho)*vbar^2*tbar^2/(2*L*t0), vbar*tbar/t0, 0 ];
% 
% A{2} = [ -rho*vbar*tbar/(L*t0),      0,        0;
%         rho*vbar*tbar/(L*t0),      0,        0;
%         rho*kappa*vbar^2*tbar^2/(2*L*t0), kappa*vbar*tbar/t0, 0 ];
% 
% Ad{2} = [ -(1-rho)*vbar*tbar/(L*t0),      0,        0;
%          (1-rho)*vbar*tbar/(L*t0),      0,        0;
%          (1-rho)*kappa*vbar^2*tbar^2/(2*L*t0), kappa*vbar*tbar/t0, 0 ];
% 
% B{1} = [ vbar*tbar/(l*t0);
%        0;
%        0 ];
% 
% B{2} = [ vbar*tbar/(l*t0);
%        0;
%        0 ];


%%%%%%%%%%%%%%%%%%% 例子2 
%Inverted pendulum example

% Parameters
M     = 1.378;
m     = 0.051;
g     = 9.8;
l     = 0.325;
cr    = 0.051;
gr    = 0.7;

alpha = cos(pi/6);
beta  = 3/pi;

% System matrices
A{1} = [ 0,      0,              1,       0;
       0,      0,              0,       1;
       0, -m*g/M,          -cr/M,       0;
       0, (M+m)*g/(M*l),   cr/(M*l),    0 ];

B{1} = [ 0;
       0;
       1/M;
      -1/(M*l) ];

A{2} = [ 0,      0,                    1,        0;
       0,      0,                    0,        1;
       0, -m*g*beta/(M*alpha),    -cr/M,       0;
       0, (M+m)*g*beta/(M*l*alpha^2), cr/(M*l*alpha), 0 ];

B{2} = [ 0;
       0;
       1/M;
      -1/(M*l*alpha) ];

Ad{1} = [ 0, 0,       0,      0;
        0, 0,       0,      0;
        0, 0,   -gr/M,      0;
        0, 0,    gr/(M*l),  0 ];

Ad{2} = [ 0, 0,       0,      0;
        0, 0,       0,      0;
        0, 0,   -gr/M,      0;
        0, 0,    gr/(M*l*alpha), 0 ];



%% system matrices (example, replace with yours)
  
[n,~] = size(A{1});
epsilon  = 1e-7;


mu1 = -mu;
mu2 = mu;


%% decision variables
P0 = sdpvar(7*n,7*n,'symmetric');
P1 = sdpvar(7*n,7*n,'symmetric');

Q1 = sdpvar(5*n,5*n,'symmetric');
Q2 = sdpvar(5*n,5*n,'symmetric');

R11 = sdpvar(n,n,'symmetric');
R12 = sdpvar(n,n,'symmetric');
R21 = sdpvar(n,n,'symmetric');
R22 = sdpvar(n,n,'symmetric');
R31 = sdpvar(n,n,'symmetric');
R32 = sdpvar(n,n,'symmetric');

Y1til = cell(1,2);
Y1til{1} = sdpvar(3*n,3*n,'full');
Y1til{2} = sdpvar(3*n,3*n,'full');

Y2til = cell(1,2);
Y2til{1} = sdpvar(3*n,3*n,'full');
Y2til{2} = sdpvar(3*n,3*n,'full');

M1til = cell(1,2);

M1til{1} = sdpvar(3*n,3*n,'symmetric');
M1til{2} = sdpvar(3*n,3*n,'symmetric');

M2til = cell(1,2);

M2til{1} = sdpvar(3*n,3*n,'symmetric');
M2til{2} = sdpvar(3*n,3*n,'symmetric');


J = cell(1,rule);
for j = 1 : rule
    J{j} = sdpvar(1,n);
end

X = sdpvar(n,n,'full');

con1 = 0.5;


%%%% line integral：
L0  = sdpvar(n,n,'symmetric');


% Li = L0 + L1_i

L = cell(1,rule);
L1 = cell(1,rule);
for i = 1:rule
    d = sdpvar(n,1);      % diagonal variables
    L1{i} = diag(d);      % P1_i is diagonal
    L{i} = L0 + L1{i};    % Pi
end



v = cell(1,10);
for i=1:10
    v{i} = [zeros(n,(i-1)*n) eye(n) zeros(n,(10-i)*n)];
end


R11T = blkdiag(R11,3*R11,5*R11);
R21T = blkdiag(R21,3*R21,5*R21);
R12T = blkdiag(R12,3*R12,5*R12);
R32T = blkdiag(R32,3*R32,5*R32);



%%%%%%%%%%%%%%% control design
%%%% 中间过程
EV0_mu1 = cell(rule, rule);
EVh2_mu1 = cell(rule, rule);
EVh_mu1 = cell(rule, rule);
EV0_mu2 = cell(rule, rule);
EVh2_mu2 = cell(rule, rule);
EVh_mu2 = cell(rule, rule);

F0_mu1 = cell(rule, rule);
F1_mu1 = cell(rule, rule);
F2_mu1 = cell(rule, rule);
F0_mu2 = cell(rule, rule);
F1_mu2 = cell(rule, rule);
F2_mu2 = cell(rule, rule);


for i  = 1 : rule
    for j = 1 : rule
EV0_mu1{i}{j}  = LMIfunc(A{i},Ad{i},B{i},0,   mu1,h,mu1,mu2,P0,P1,Q1,Q2,R11,R12,R21,R22,R31,R32,L{i},Y1til,Y2til,M1til,M2til,con1,X,J{j},n);
EVh2_mu1{i}{j} = LMIfunc(A{i},Ad{i},B{i},h/2, mu1,h,mu1,mu2,P0,P1,Q1,Q2,R11,R12,R21,R22,R31,R32,L{i},Y1til,Y2til,M1til,M2til,con1,X,J{j},n);
EVh_mu1{i}{j}  = LMIfunc(A{i},Ad{i},B{i},h,   mu1,h,mu1,mu2,P0,P1,Q1,Q2,R11,R12,R21,R22,R31,R32,L{i},Y1til,Y2til,M1til,M2til,con1,X,J{j},n);

EV0_mu2{i}{j}  = LMIfunc(A{i},Ad{i},B{i},0,   mu2,h,mu1,mu2,P0,P1,Q1,Q2,R11,R12,R21,R22,R31,R32,L{i},Y1til,Y2til,M1til,M2til,con1,X,J{j},n);
EVh2_mu2{i}{j} = LMIfunc(A{i},Ad{i},B{i},h/2, mu2,h,mu1,mu2,P0,P1,Q1,Q2,R11,R12,R21,R22,R31,R32,L{i},Y1til,Y2til,M1til,M2til,con1,X,J{j},n);
EVh_mu2{i}{j}  = LMIfunc(A{i},Ad{i},B{i},h,   mu2,h,mu1,mu2,P0,P1,Q1,Q2,R11,R12,R21,R22,R31,R32,L{i},Y1til,Y2til,M1til,M2til,con1,X,J{j},n);

    end
end

for i = 1  : rule
    for j = 1 : rule
    F0_mu1{i}{j} = EV0_mu1{i}{j};
    F1_mu1{i}{j} = (4*EVh2_mu1{i}{j} - EVh_mu1{i}{j} - 3*EV0_mu1{i}{j})/h;
    F2_mu1{i}{j} = 2*(EVh_mu1{i}{j} - 2*EVh2_mu1{i}{j} + EV0_mu1{i}{j})/h^2;

    F0_mu2{i}{j} = EV0_mu2{i}{j};
    F1_mu2{i}{j} = (4*EVh2_mu2{i}{j} - EVh_mu2{i}{j} - 3*EV0_mu2{i}{j})/h;
    F2_mu2{i}{j}= 2*(EVh_mu2{i}{j} - 2*EVh2_mu2{i}{j} + EV0_mu2{i}{j})/h^2;
    end
end



%%%% 比较对象
Sur11 = [R11T + R21T  , -Y1til{1};
        -Y1til{1}'  ,  M1til{1}];

Sur12 = [R11T + R21T  , -Y1til{2};
        -Y1til{2}'  ,  M1til{2}];

Sur21 = [M2til{1}      , -Y2til{1};
        -Y2til{1}'  , R12T + R32T];

Sur22 = [M2til{2}      , -Y2til{2};
        -Y2til{2}' , R12T + R32T];

LMI1 = Sur11;   
LMI2 = Sur12;
LMI3 = Sur21;   
LMI4 = Sur22;

LMIcon1 = cell(1,rule);
LMIcon2 = cell(1,rule);
LMIcon3 = cell(1,rule);
LMIcon4 = cell(1,rule);
LMIcon5 = cell(1,rule);
LMIcon6 = cell(1,rule);
LMIcon7 = cell(1,rule);
LMIcon8 = cell(1,rule);

c1 = cell(1,rule);
c1{1} = 0.7;
c1{2} = 0.7;


c2 = cell(1,rule);
c2{1}{2} = 0.7; 

%%%% i = j
for i = 1 :rule
    LMIcon1{i} = F0_mu1{i}{i};
    LMIcon2{i} = F0_mu2{i}{i};
 
    LMIcon3{i} = h^2* F2_mu1{i}{i} + h* F1_mu1{i}{i} + F0_mu1{i}{i};
    LMIcon4{i} = h^2* F2_mu2{i}{i} + h* F1_mu2{i}{i} + F0_mu2{i}{i};

    LMIcon5{i} = c1{i}*h* F1_mu1{i}{i} + 2* F0_mu1{i}{i};
    LMIcon6{i} = c1{i}*h* F1_mu2{i}{i}  + 2* F0_mu2{i}{i};

    LMIcon7{i} = 2*c1{i}*h^2* F2_mu1{i}{i} + (1 + c1{i})*h* F1_mu1{i}{i} + 2*F0_mu1{i}{i};
    LMIcon8{i} = 2*c1{i}*h^2* F2_mu2{i}{i} + (1 + c1{i})*h* F1_mu2{i}{i} + 2*F0_mu2{i}{i};
end



LMI2con1 = cell(rule,rule-1);
LMI2con2 = cell(1,rule);
LMI2con3 = cell(1,rule);
LMI2con4 = cell(1,rule);
LMI2con5 = cell(1,rule);
LMI2con6 = cell(1,rule);
LMI2con7 = cell(1,rule);
LMI2con8 = cell(1,rule);


%%%% j > i 
for i = 1 :rule
    for j = i + 1 : rule  
    LMI2con1{i}{j} = F0_mu1{i}{j} + F0_mu1{j}{i};
    LMI2con2{i}{j} = F0_mu2{i}{j} + F0_mu2{j}{i};

    LMI2con3{i}{j} = h^2* (F2_mu1{i}{j} + F2_mu1{j}{i}) + h* (F1_mu1{i}{j} + F1_mu1{j}{i})  + (F0_mu1{i}{j}+F0_mu1{j}{i});
    LMI2con4{i}{j} = h^2* (F2_mu2{i}{j} + F2_mu2{j}{i}) + h* (F1_mu2{i}{j} + F1_mu2{j}{i})  + (F0_mu2{i}{j}+F0_mu2{j}{i});

    LMI2con5{i}{j} = c2{i}{j}*h* (F1_mu1{i}{j} + F1_mu1{j}{i})  + 2* (F0_mu1{i}{j}+F0_mu1{j}{i});
    LMI2con6{i}{j} = c2{i}{j}*h* (F1_mu2{i}{j} + F1_mu2{j}{i})  + 2* (F0_mu2{i}{j}+F0_mu2{j}{i});

    LMI2con7{i}{j} = 2*c2{i}{j}*h^2* (F2_mu1{i}{j} + F2_mu1{j}{i}) + (1 + c2{i}{j})*h* (F1_mu1{i}{j} + F1_mu1{j}{i}) + 2*(F0_mu1{i}{j}+F0_mu1{j}{i});
    LMI2con8{i}{j} = 2*c2{i}{j}*h^2* (F2_mu2{i}{j} + F2_mu2{j}{i}) + (1 + c2{i}{j})*h* (F1_mu2{i}{j} + F1_mu2{j}{i}) + 2*(F0_mu2{i}{j}+F0_mu2{j}{i});
    end
end





constraints = [
        LMI1 >= 0;
        LMI2 >= 0;
        LMI3 >= 0;
        LMI4 >= 0;
        P0 >= epsilon * eye(7*n);
        P0+ h*P1 >= epsilon * eye(7*n);
        Q1 >= epsilon * eye(5*n);
        Q2 >= epsilon * eye(5*n);
        R11 >= epsilon * eye(n);
        R12 >= epsilon * eye(n);
        R21 >= epsilon * eye(n);
        R22 >= epsilon * eye(n);
        R31 >= epsilon * eye(n);
        R32 >= epsilon * eye(n);
                  ];

constraints = [constraints, diag(L0) == 0 ];

for i = 1 : rule
    constraints = [constraints, LMIcon1{i}<= -epsilon*eye(10*n), 
                                LMIcon2{i}<= -epsilon*eye(10*n), 
                                LMIcon3{i}<= -epsilon*eye(10*n), 
                                LMIcon4{i}<= -epsilon*eye(10*n), 
                                LMIcon5{i}<= -epsilon*eye(10*n), 
                                LMIcon6{i}<= -epsilon*eye(10*n), 
                                LMIcon7{i}<= -epsilon*eye(10*n), 
                                LMIcon8{i}<= -epsilon*eye(10*n), 
                                 L{i} >= epsilon*eye(n), ];
end


for i = 1 : rule
    for j = i + 1 : rule 
    constraints = [constraints, LMI2con1{i}{j}<= -epsilon*eye(10*n), 
                                LMI2con2{i}{j}<= -epsilon*eye(10*n), 
                                LMI2con3{i}{j}<= -epsilon*eye(10*n), 
                                LMI2con4{i}{j}<= -epsilon*eye(10*n), 
                                LMI2con5{i}{j}<= -epsilon*eye(10*n), 
                                LMI2con6{i}{j}<= -epsilon*eye(10*n), 
                                LMI2con7{i}{j}<= -epsilon*eye(10*n), 
                                LMI2con8{i}{j}<= -epsilon*eye(10*n), ];
    end
end



    options = sdpsettings('solver', 'mosek', 'verbose', 0);
    % options = sdpsettings('solver', 'sdpt3', 'verbose', 0);




solll = optimize(constraints, [], options);

if solll.problem == 0

    feasibli = 0;

    X_value = value(X);

    J1 = value(J{1});

    J2 = value(J{2});

    K1 = J1/X_value;

    K2 = J2/X_value;

    fprintf('X =\n');

    disp(X_value);

    fprintf('J1 =\n');

    disp(J1);

    fprintf('J2 =\n');

    disp(J2);

    fprintf('K1 =\n');

    disp(K1);

    fprintf('K2 =\n');

    disp(K2);

else

    feasibli = 1;

end
end



function [EE]=LMIfunc(A,Ad,B,pit,pit_dot,h,mu1,mu2,P0,P1,Q1,Q2,R11,R12,R21,R22,R31,R32,L,Y1til,Y2til,M1til,M2til,con1,X,J,n)

v = cell(1,10);
for i=1:10
    v{i} = [zeros(n,(i-1)*n) eye(n) zeros(n,(10-i)*n)];
end

%% auxiliary vectors
v0 = zeros(n,10*n);
alpha = pit/h;
m21 = mu2 - mu1;
Y1dot =  (mu2 - pit_dot)/m21 * Y1til{1} + (pit_dot - mu1)/m21 * Y1til{2};
Y2dot =  (mu2 - pit_dot)/m21 * Y2til{1} + (pit_dot - mu1)/m21 * Y2til{2};

M1dot =  (mu2 - pit_dot)/m21 * M1til{1} + (pit_dot - mu1)/m21 * M1til{2};
M2dot =  (mu2 - pit_dot)/m21 * M2til{1} + (pit_dot - mu1)/m21 * M2til{2};

%% G1
G1 = [v{1}',v{2}',v{3}',(pit*v{6})',((h-pit)*v{7})',(pit*v{8})',((h-pit)*v{9})']';

%% G2
G2 = [v{10}',((1-pit_dot)*v{4})',v{5}',...
         (v{1}-(1-pit_dot)*v{2})',...
         ((1-pit_dot)*v{2}-v{3})',...
         (v{1}-(1-pit_dot)*v{6}-pit_dot*v{8})',...
         ((1-pit_dot)*v{2}-v{7}+pit_dot*v{9})']';

%% G3
G3 = [v{10}',v{1}',v{1}',v0',(pit*v{6})']';

%% G4
G4 = [v{4}',v{2}',v{1}',(pit*v{6})',v0']';

%% G5
G5 = [v{4}',v{2}',v{1}',v0',((h-pit)*v{7})']';

%% G6
G6 = [v{5}',v{3}',v{1}',((h-pit)*v{7})',v0']';

%% G7
G7 = [(v{1}-v{2})',(pit*v{6})',(pit*v{1})',(pit^2*v{8})',(pit^2*v{6}-pit^2*v{8})']';

%% G8
G8 = [v0',v0', v{10}',v{1}',(-(1-pit_dot)*v{2})']';

%% G9
G9 = [(v{2}-v{3})',((h-pit)*v{7})', ((h-pit)*v{1})',((h-pit)^2*v{9})',...
         ((h-pit)^2*v{7}-(h-pit)^2*v{9})']';

%% G10
G10 = [v0',v0',v{10}',((1-pit_dot)*v{2})',(-v{3})']';

%% G11
G11 = [v{1}' - v{6}', v{1}' + 2*v{6}' - 6*v{8}' ]';
%% G12
G12 = [v{2}' - v{6}', v{2}' - 4*v{6}' + 6*v{8}' ]';
%% G13
G13 = [v{2}' - v{7}', v{2}' + 2*v{7}' - 6*v{9}' ]';
%% G14
G14 = [v{3}' - v{7}', v{3}' - 4*v{7}' + 6*v{9}' ]';

%% Omega1 ; Omega2
Omega1 = [(v{1}-v{2})',(v{1}+v{2}-2*v{6})',(v{1}-v{2}+6*v{6}-12*v{8})']';
Omega2 = [(v{2}-v{3})',(v{2}+v{3}-2*v{7})',(v{2}-v{3}+6*v{7}-12*v{9})']';

%% Rm
R11T = blkdiag(R11,3*R11,5*R11);
R12T = blkdiag(R12,3*R12,5*R12);
R21T = blkdiag(R21,3*R21,5*R21);
R32T = blkdiag(R32,3*R32,5*R32);


R21S = blkdiag(2*R21,4*R21);
R22S = blkdiag(2*R22,4*R22);
R31S = blkdiag(2*R31,4*R31);
R32S = blkdiag(2*R32,4*R32);


%% PN
PN = P0 + pit*P1;
dPN =  pit_dot*P1;

%% Π1
Pi1 = He(G1'*PN*G2) + G1'*(dPN)*G1;

%% Π2
Pi2 = G3'*Q1*G3 ...
     -(1-pit_dot)*G4'*Q1*G4 ...
     +(1-pit_dot)*G5'*Q2*G5 ...
     -G6'*Q2*G6 ...
     +He(G7'*Q1*G8 + G9'*Q2*G10);

%% Π3

m = 3;
sum1 = 0;
sum2 = 0;
for i = 1 : m-1
    sum1 = sum1 + (1-alpha)^(i)*(R11T+R21T);
    sum2 = sum2 +  alpha^(i)*(R12T+R32T);
end


phi1 =  R11T + sum1;
phi2 =  (1-alpha)^(m-1)*Y1dot + alpha^(m-1)*Y2dot ;
phi3 =  R12T + sum2;

phv1 =  alpha^(m-2)*(1-alpha)*M2dot;
phv2 =  (alpha)*(1-alpha)^(m-2)*M1dot;




Pi3 = h^2*v{10}'*(2*R11 + R21 + R31)*v{10} ...
      - (1 - pit_dot)*(h - pit)* v{4}'*(2*h*(R11 - R12) + (h - pit)*(R21 - R22) + (h + pit)*(R31 - R32))*v{4};
%% Π4


%% 
%% 
Pi4 = -2*(G11'*R21S*G11+G12'*R31S*G12+ G13'*R22S*G13 + G14'*R32S*G14 )  ...
      -2*[Omega1;Omega2]'*[phi1-phv1, phi2; 
                           phi2', phi3- phv2]*[Omega1;Omega2];

%% Π5 控制器设计：

Pi5 = (v{1}' + con1*v{10}') * (-X*v{10} + (A*X + B*J) * v{1} + Ad*X*v{2});

%% Π6(line-integral)
Pi6 = v{1}'* L * v{10}; 

%% total LMI
EE = Pi1 + Pi2 + Pi3 + Pi4 + He(Pi5 + Pi6);

end


function SYM = He(X)
SYM = X+X';
end