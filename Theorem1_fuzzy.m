clear
clc
% 
yalmip('clear')
mu_values = [0.2, 0.4, 0.6];  %% example one
%mu_values = [1e-7,0.1, 0.5];   %% example two
tol = 1e-4;                        

%%%% bisection method
    for mu_idx = 1:length(mu_values)
        mu = mu_values(mu_idx);
        fprintf('Processing mu=%.1f\n', mu);
        pi_min = 0;
        pi_max = 18;
        while pi_max - pi_min > tol
            h = (pi_min + pi_max) / 2;
            feasibli = check_lmi(mu, h);
            if feasibli == 0 
                pi_min = h;  
            else
                pi_max = h;  
            end
        end
        fprintf('mu=%.1f, tau=%.4f\n', mu, pi_min);
    end



function feasibli = check_lmi(mu,h)

%% system dimension

rule = 2;

A = cell(1, rule);
Ad = cell(1, rule);


%% example one
A{1} = [-2   0; 
        0 -0.9];
Ad{1} = [-1 0;
        -1 -1];

A{2} = [-1.5  1;
         0   -0.75];
Ad{2} = [-1   0;
          1  -0.85];


% %%%%% example two
% A{1} = [-2.1   0.1; 
%         -0.2   -0.9];
% 
% Ad{1} = [-1.1  0.1;
%         -0.8   -0.9];
% 
% A{2} = [-1.9   0;
%         -0.2   -1.1];
% 
% Ad{2} = [-0.9    0;
%          -1.1   -1.2];


%% system matrices (example)
  
[n,~] = size(A{1});
epsilon  = 1e-7; %%% tolerance of Yalmip

%%% example 
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


K1 = sdpvar(n,n,'full');
K2 = sdpvar(n,n,'full');
K3 = sdpvar(n,n,'full');


%%%%%%% line integral：
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

%%%%  the process to calculate coefficient of quadratic polynomial
EV0_mu1 = cell(1, rule);
EVh2_mu1 = cell(1, rule);
EVh_mu1 = cell(1, rule);
EV0_mu2 = cell(1, rule);
EVh2_mu2 = cell(1, rule);
EVh_mu2 = cell(1, rule);

F0_mu1 = cell(1, rule);
F1_mu1 = cell(1, rule);
F2_mu1 = cell(1, rule);
F0_mu2 = cell(1, rule);
F1_mu2 = cell(1, rule);
F2_mu2 = cell(1, rule);


for i  = 1 : rule
EV0_mu1{i}  = LMIfunc(A{i},Ad{i},0,   mu1,h,mu1,mu2,P0,P1,Q1,Q2,R11,R12,R21,R22,R31,R32,L{i},Y1til,Y2til,M1til,M2til,K1,K2,K3,n);
EVh2_mu1{i} = LMIfunc(A{i},Ad{i},h/2, mu1,h,mu1,mu2,P0,P1,Q1,Q2,R11,R12,R21,R22,R31,R32,L{i},Y1til,Y2til,M1til,M2til,K1,K2,K3,n);
EVh_mu1{i}  = LMIfunc(A{i},Ad{i},h,   mu1,h,mu1,mu2,P0,P1,Q1,Q2,R11,R12,R21,R22,R31,R32,L{i},Y1til,Y2til,M1til,M2til,K1,K2,K3,n);

EV0_mu2{i}  = LMIfunc(A{i},Ad{i},0,   mu2,h,mu1,mu2,P0,P1,Q1,Q2,R11,R12,R21,R22,R31,R32,L{i},Y1til,Y2til,M1til,M2til,K1,K2,K3,n);
EVh2_mu2{i} = LMIfunc(A{i},Ad{i},h/2, mu2,h,mu1,mu2,P0,P1,Q1,Q2,R11,R12,R21,R22,R31,R32,L{i},Y1til,Y2til,M1til,M2til,K1,K2,K3,n);
EVh_mu2{i}  = LMIfunc(A{i},Ad{i},h,   mu2,h,mu1,mu2,P0,P1,Q1,Q2,R11,R12,R21,R22,R31,R32,L{i},Y1til,Y2til,M1til,M2til,K1,K2,K3,n);
end




for i = 1  : rule
    F0_mu1{i} = EV0_mu1{i};
    F1_mu1{i} = (4*EVh2_mu1{i} - EVh_mu1{i} - 3*EV0_mu1{i})/h;
    F2_mu1{i} = 2*(EVh_mu1{i} - 2*EVh2_mu1{i} + EV0_mu1{i})/h^2;

    F0_mu2{i} = EV0_mu2{i};
    F1_mu2{i} = (4*EVh2_mu2{i} - EVh_mu2{i} - 3*EV0_mu2{i})/h;
    F2_mu2{i}= 2*(EVh_mu2{i} - 2*EVh2_mu2{i} + EV0_mu2{i})/h^2;
end



%%%% write LMIs
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


%%%%%%%%%%%%%% value c^i please change number
c = cell(1,rule);
c{1} = 0.7;
c{2} = 0.7;



for i = 1 :rule
    LMIcon1{i} = F0_mu1{i};
    LMIcon2{i} = F0_mu2{i};

    LMIcon3{i} = h^2* F2_mu1{i} + h* F1_mu1{i} + F0_mu1{i};
    LMIcon4{i} = h^2* F2_mu2{i} + h* F1_mu2{i} + F0_mu2{i};

    LMIcon5{i} = c{i}*h* F1_mu1{i}  + 2* F0_mu1{i};
    LMIcon6{i} = c{i}*h* F1_mu2{i}  + 2* F0_mu2{i};

    LMIcon7{i} = 2*c{i}*h^2* F2_mu1{i} + (1 + c{i})*h* F1_mu1{i} + 2*F0_mu1{i};
    LMIcon8{i} = 2*c{i}*h^2* F2_mu2{i} + (1 + c{i})*h* F1_mu2{i} + 2*F0_mu2{i};
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

%%% line integral constraint
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
                                 L{i} >= epsilon*eye(n)];
end

    options = sdpsettings('solver', 'mosek', 'verbose', 0);
   




solll = optimize(constraints, [], options);

if solll.problem == 0
    feasibli = 0;
else
    feasibli = 1;

end
end



function [EE]=LMIfunc(A,Ad,pit,pit_dot,h,mu1,mu2,P0,P1,Q1,Q2,R11,R12,R21,R22,R31,R32,L,Y1til,Y2til,M1til,M2til,W1,W2,W3,n)

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
X1 = [(v{1}-v{2})',(v{1}+v{2}-2*v{6})',(v{1}-v{2}+6*v{6}-12*v{8})']';
X2 = [(v{2}-v{3})',(v{2}+v{3}-2*v{7})',(v{2}-v{3}+6*v{7}-12*v{9})']';

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

%%%% please change
%m = 2;
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


Pi4 = -2*(G11'*R21S*G11+G12'*R31S*G12+ G13'*R22S*G13 + G14'*R32S*G14 )  ...
      -2*[X1;X2]'*[phi1-phv1, phi2; 
                           phi2', phi3- phv2]*[X1;X2];


%%%% line integral
Pi5 = v{1}'* L * v{10} ; 
%% Π6 zero equation
Pi6 = (v{1}'*W1 + v{2}'*W2 + v{10}'*W3)* (-v{10} + A * v{1} + Ad* v{2});


%% total LMI
EE = Pi1 + Pi2 + Pi3 + Pi4 + He(Pi5 + Pi6);

end


function SYM = He(X)
SYM = X+X';
end
