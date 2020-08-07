% Main code for Van-Dinh Nguyen, Le-Nam Tran, Trung Q. Duong, Oh-Soon Shin, and Ronan Farrell,
%"An Efficient Precoder Design for Multiuser MIMO Cognitive Radio Networks with Interference Constraints," 
%IEEE Transactions on Vehicular Technology, vol. 66, no. 5, pp. 3991-4004, May 2017.



clear all
clc
nTx = 10; % number of users
nRxarray = [2 2 2 2 2];   %  number of receive antennas for each user, i.e.,
                        %  nRxarray(i) is the number of rx. antennas at the ith user
P_dB = [0:2:20]; % sum power, dB scale                    
nUsers = length(nRxarray); % number of users
finalsumrate=zeros(length(P_dB),1);
%%%%%%%%%%%%%%%%%

nRxSU = [2 2 2]; % nRxSU(k) is the # of rx antenna at SU k.
nRXPU = [2 2 2 2 2]; % nRxPU(k) is the # of rx antenna at PU k.
totalpowermax_dB =10; %max total power in dB
Interferencetheshold_dB=1;
totalpowermax = 10^(totalpowermax_dB/10); 
nSU = length(nRxSU); % # of SUs
nPU = length(nRXPU); % # of PUs
S = cell(nSU,1);
Interferencetheshold = (10^(Interferencetheshold_dB/10))*ones(nPU,1);
channelSU = cell(nSU,1);
effectivechannelSU = cell(nSU,1);
channelPU = cell(nPU,1);
nullspacebasis = cell(nSU,1);

%%%%%%%%%%%%%%%%%%%%


for nn=1:length(P_dB)
P = 10.^(P_dB(nn)./10);
Po = P/nTx*ones(nTx,1); % maximum power per antenna 
Numberrunning = 1000;% The number of running times
sumraterandom=zeros(Numberrunning,1);
channel = cell(nUsers,1); % channel matrices of users
for n=1:Numberrunning
nn
n
powerperantmax = totalpowermax/nTx; % max  power per antenna
for iSU = 1:nSU
    channelSU{iSU} = 1/sqrt(2)*(randn(nRxSU(iSU),nTx) + 1i*randn(nRxSU(iSU),nTx));
end
for iPU=1:nPU
    channelPU{iPU}  = 1/sqrt(2)*(randn(nRXPU(iPU),nTx) + 1i*randn(nRXPU(iPU),nTx));
end


F = [];
sumrateSU = 0;
powerperant = 0;
for iSU = 1:nSU
    tempchan = cell2mat(channelSU((1:nSU~=iSU)));
    nullspacebasis{iSU}= null(tempchan);
    effectivechannelSU{iSU} = channelSU{iSU}*nullspacebasis{iSU};
   % S{iSU} = sdpvar(size(nullspacebasis{iSU},2),size(nullspacebasis{iSU},2),'hermitian','complT = diag(mylambda);ex');
   S{iSU} = sdpvar(size(nullspacebasis{iSU},2),size(nullspacebasis{iSU},2),'hermitian','complex');
    F = [F,S{iSU} >= 0];
    sumrateSU = sumrateSU + logdet(eye(nRxSU(iSU))+effectivechannelSU{iSU}*S{iSU}*...
        effectivechannelSU{iSU}');
    
    powerperant = powerperant + nullspacebasis{iSU}*S{iSU}*(nullspacebasis{iSU})';
end
F = [F,real(diag(powerperant)) <= powerperantmax];

for iPU = 1:nPU
    F = [F, real(channelPU{iPU}*powerperant*channelPU{iPU}') <= eye(2)*Interferencetheshold(iPU)];
    %F = [F, real(min(eig(channelPU{iPU}*powerperant*channelPU{iPU}'))) >= 0];
end
  %  F = [F,real(trace(channelPU{nPU}*powerperant*channelPU{nPU}')) <= 10];
myobj = sumrateSU;
myops = sdpsettings('solver','sedumi','verbose',0);
diagnotics = solvesdp(F,-myobj,myops);
diagnotics.solvertime;
temp=double(myobj);
precodersoftwareSU = cell(nSU,1);
 S_sum=0;
for iSU = 1:nSU
     precodersoftwareSU{iSU} = double(S{iSU});
    S_sum=S_sum+precodersoftwareSU{iSU};
end
INT=zeros(nPU,1);
for iPU = 1:nPU
  INT(iPU)= trace(channelPU{iPU}*nullspacebasis{iSU}*S_sum*nullspacebasis{iSU}'*channelPU{iPU}');
end



for iUser=1:nUsers
    channel{iUser} = sqrt(1/2)*(randn(nRxarray(iUser),nTx) + 1i*randn(nRxarray(iUser),nTx));    
end

%semilogy(gap)
%% YALMIP code for solving (31)
effectivechannel = cell(nUsers,1);
effectivechannel{1} = channel{1};
B = cell(nUsers,1);
B{1} = eye(nTx);
h_bar = [];
A = cell(nUsers,nTx); % matrix A in (33)
for iUser=1:nUsers
    h_bar = [h_bar;channel{iUser}];
    if( iUser < nUsers)
        B{iUser+1} = null(h_bar);
    end
    effectivechannel{iUser} = channel{iUser}*B{iUser};
    for iTx=1:nTx
        A{iUser,iTx}=B{iUser}'*diag(1:nTx==iTx)*B{iUser}; % compute matrix A in (33)
    end    
end


X=cell(nUsers,1);
F=[];
for iUser = 1:nUsers; 
    X{iUser} = sdpvar(nTx-sum(nRxarray(1:iUser-1)),nTx-sum(nRxarray(1:iUser-1)),'hermitian','complex');
    F=[F,X{iUser}>=0];
end
obj=0;
for iUser=1:nUsers
    obj=obj+logdet(eye(nRxarray(iUser)) + (effectivechannel{iUser}*X{iUser}*effectivechannel{iUser}')/(1+INT(iUser)));
end

for iTx=1:nTx
    y = 0;
    for iUser=1:nUsers
        y = y + real(trace(A{iUser,iTx}*X{iUser}));
    end
    F=[F,y <= Po(iTx)];
end

ops = sdpsettings('solver','sdpt3','verbose',0);

sol = solvesdp(F,-obj,ops);

sumrateyalmip = real(double(obj));
precodersoftware = cell(nUsers,1);
for iUser = 1:nUsers
    precodersoftware{iUser} = B{iUser}*double(X{iUser})*B{iUser}';
end
sumraterandom(n)=sumrateyalmip;
end
finalsumrate(nn)=real(mean(sumraterandom));
end
%plot(totalpowermax_dB, finalsumrate);
plot(P_dB, finalsumrate, 'k^--', 'linewidth', 2, 'markersize',7);
