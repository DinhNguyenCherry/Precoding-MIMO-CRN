% Main code for Van-Dinh Nguyen, Le-Nam Tran, Trung Q. Duong, Oh-Soon Shin, and Ronan Farrell,
%"An Efficient Precoder Design for Multiuser MIMO Cognitive Radio Networks with Interference Constraints," 
%IEEE Transactions on Vehicular Technology, vol. 66, no. 5, pp. 3991-4004, May 2017.


clear all
clc
nTx  = 10;%number of antennas at the BS-SU
nRxSU = [2 2 2]; % nRxSU(k) is the # of rx antenna at SU k.
nRXPU = [2 2]; % nRxPU(k) is the # of rx antenna at PU k.
totalpowermax_dB =15:3:15; %max total power in dB
finalsumrate=zeros(length(totalpowermax_dB),1);

for nn=1:length(totalpowermax_dB)%running in ranging of transmitpower
    totalpowermax = 10^(totalpowermax_dB(nn)/10);%convert from dB to power  
    %powerperantmax = totalpowermax/nTx; % max  power per antenna
    nSU = length(nRxSU); % # of SUs
    nPU = length(nRXPU); % # of PUs
    S = cell(nSU,1);
    Interferencetheshold = (10^(10/10))*ones(nPU,1);%transmit power constraint
    channelSU = cell(nSU,1);
    effectivechannelSU = cell(nSU,1);
    channelPU = cell(nPU,1);
    nullspacebasis = cell(nSU,1);
    Numberrunning = 1;% The number of running times to compute average sum rate
    sumraterandom=zeros(Numberrunning,1);%average sum rate over Numberofrunning

   for n=1:Numberrunning
 %% The optimization for problem of (4)
           nn
           n
           powerperantmax = totalpowermax/nTx; % max  power per antenna
        for iSU = 1:nSU % Generate the channel of SU
            channelSU{iSU} = 1/sqrt(2)*(randn(nRxSU(iSU),nTx) + 1i*randn(nRxSU(iSU),nTx));
        end
        for iPU=1:nPU % Generate the channel of PU
            channelPU{iPU}  = 1/sqrt(2)*(randn(nRXPU(iPU),nTx) + 1i*randn(nRXPU(iPU),nTx));
        end
        sumrate = 0;
        F = [];
        sumrate = 0;
        powerperant = 0;
        for iSU = 1:nSU
            tempchan = cell2mat(channelSU((1:nSU~=iSU)));
            nullspacebasis{iSU}= null(tempchan);
            effectivechannelSU{iSU} = channelSU{iSU}*nullspacebasis{iSU};
            S{iSU} = sdpvar(size(nullspacebasis{iSU},2),size(nullspacebasis{iSU},2),'hermitian','complex');
            F = [F,S{iSU} >= 0];
            sumrate = sumrate + logdet(eye(nRxSU(iSU))+effectivechannelSU{iSU}*S{iSU}*...
            effectivechannelSU{iSU}');
            powerperant = powerperant + nullspacebasis{iSU}*S{iSU}*(nullspacebasis{iSU})';
        end
         F = [F,real(diag(powerperant)) <= powerperantmax];
        for iPU = 1:nPU
            F = [F,real(trace(channelPU{iPU}*powerperant*channelPU{iPU}')) <= Interferencetheshold(iPU)];
        end
    obj = sumrate;
    myops = sdpsettings('solver','sedumi','verbose',0);
    diagnotics = solvesdp(F,-obj,myops);
    diagnotics.solvertime;
    A1=double(obj);% sum rate of (4)


%% Transfer to real domain
    Sreal = cell(nSU,1);
    sumrate = 0;
    F = [];
    powerperant = 0;
        for iSU = 1:nSU
            Sreal{iSU} = sdpvar(2*size(nullspacebasis{iSU},2),2*size(nullspacebasis{iSU},2));
            F = [F,Sreal{iSU} >= 0];
            sumrate = sumrate + 1/2*logdet(eye(2*nRxSU(iSU))+...
                [real(effectivechannelSU{iSU}) -imag(effectivechannelSU{iSU});
            imag(effectivechannelSU{iSU}) real(effectivechannelSU{iSU})]*Sreal{iSU}*...
            [real(effectivechannelSU{iSU}) -imag(effectivechannelSU{iSU});
            imag(effectivechannelSU{iSU}) real(effectivechannelSU{iSU})]');

            powerperant = powerperant + [(nullspacebasis{iSU}),1i*(nullspacebasis{iSU})]...
                *(Sreal{iSU})*[(nullspacebasis{iSU}),1i*(nullspacebasis{iSU})]';
        end
            F = [F,1/2*real((diag(powerperant))) <= powerperantmax];
        for iPU = 1:nPU
            F = [F,1/2*real(trace([(channelPU{iPU})]*powerperant*...
                [(channelPU{iPU})]')) <= Interferencetheshold(iPU)];
        end
        obj = sumrate;
        diagnotics = solvesdp(F,-obj,myops);
        diagnotics.solvertime
        % Extract dual variables
        A2=double(obj);
        mylambda = dual(F(nSU+1));
        mynu = dual(F(nSU+2:end));


%% Dual problem in (5)
            p = [powerperantmax*ones(nTx,1);Interferencetheshold];
            mynewlambda = [mylambda;mynu];
            P = p'*mynewlambda;
            T = diag(mylambda);
            %T1 = diag(mynewlambda);
            for iPU = 1:nPU
                T = T+mynu(iPU)*(channelPU{iPU})'*channelPU{iPU};
            end
            sumratedual = 0;
            Sdual = cell(nSU,1);
            totalpowerdual = 0;
            F =[];
            for iSU = 1:nSU
                Sdual{iSU} = sdpvar(size(channelSU{iSU},1),size(channelSU{iSU},1),'hermitian','complex');
                F = [F,Sdual{iSU} >= 0];
                sumratedual = sumratedual + logdet((nullspacebasis{iSU})'*T*nullspacebasis{iSU}...
                    +(effectivechannelSU{iSU})'*Sdual{iSU}*...
                    effectivechannelSU{iSU});
                totalpowerdual = totalpowerdual + real(trace(Sdual{iSU}));
            end
            F = [F,totalpowerdual <= P];
            %F = [F,p'*mynewlambda <= totalpowermax];

            obj = sumratedual;
            myops = sdpsettings('solver','sedumi','verbose',0);
            %myops = sdpsettings('solver','sedumi','verbose',0,'dualize',1);
            diagnotics = solvesdp(F,-obj,myops);
            double(obj);
            diagnotics.solvertime
            sumratedual = double(obj);
            for iSU = 1:nSU
                sumratedual = sumratedual - log(det((nullspacebasis{iSU})'*T*nullspacebasis{iSU}));
            end
            A3=sumratedual;


%% Dual problem in real domain
            sumratedual = 0;
            Sdualreal = cell(nSU,1);
            totalpowerdual = 0;
            F =[];
            for iSU = 1:nSU
                Sdual{iSU} = sdpvar(size(channelSU{iSU},1),size(channelSU{iSU},1),'hermitian','complex');
                F = [F,Sdual{iSU} >= 0];
                sumratedual = sumratedual + logdet((nullspacebasis{iSU})'*T*nullspacebasis{iSU}...
                    +(effectivechannelSU{iSU})'*Sdual{iSU}*...
                    effectivechannelSU{iSU});
                totalpowerdual = totalpowerdual + real(trace(Sdual{iSU}));
            end
            F = [F,totalpowerdual <= P];
            obj = sumratedual;
            myops = sdpsettings('solver','sedumi','verbose',0);
            %myops = sdpsettings('solver','sedumi','verbose',0,'dualize',1);
            diagnotics = solvesdp(F,-obj,myops);
            sumratedualtrue=double(obj);
            diagnotics.solvertime;
            for iSU = 1:nSU
                sumratedualtrue = sumratedualtrue - log(det((nullspacebasis{iSU})'*T*nullspacebasis{iSU}));
            end
            A4=sumratedualtrue;
            sumraterandom(n)=sumratedualtrue;
   end
         finalsumrate(nn)=mean(sumraterandom);
end
%plot(totalpowermax_dB, finalsumrate);
plot(totalpowermax_dB, finalsumrate, 'k^--', 'linewidth', 2, 'markersize',7);
    