% Ahmet Ali Tilkicioğlu - 210102002163 - Digital Communication

clc;
clear all;
close all;

% a side code 

N = 10^7;                           % 10 million bit size
Eb = 1;                             % average bit energy
Ps1 = 0.5;                          % Probability bit 1
Ps2 = 0.5;                          % Probability bit 0
si = randsrc(1,N,[1,0;Ps1,Ps2]);    % si values 
a1 = 1;                             % bit 1 signal value
a2 = -1;                            % bit 0 signal value
ai = 2*si - 1;                      % all ai signal values   
max_dB = 17;                        % db value maximum limit   


siAllValue = zeros(1,N);            % si vector
simPbValue = zeros(0,max_dB);       % simulation vector
AnltyPbValue = zeros(0,max_dB);     % analytic summaries vector



for db_SNR = 0:max_dB               % loop for all process (0:18 db)
    
   SNR = 10^(db_SNR/10);            % SNR values 
   N0 = Eb/SNR;                     % No values
   sigma0 = sqrt(N0);               % sigma value

   z = ai + sigma0*randn(1,N);      % z = (ai + no) equations 
   gamma0 = 0;                      % gamma formula
    
   siAllValue = z > gamma0;         % decision circuit simulation
   
   errorVal = si ~= siAllValue;     % Total error for simulation 
   Pb = sum(errorVal)/N;            % Pb for simulation
   simPbValue = [simPbValue Pb];    % Pb simulation vector for graph

   analytical_sol = qfunc(1/sigma0);    % analytic solution 
   AnltyPbValue = [AnltyPbValue analytical_sol]; % analytic solution 

end

% Solution graphs code for a)

dB_SNR = 0:max_dB;
figure()
semilogy(dB_SNR,simPbValue,'go-');
hold on;
semilogy(dB_SNR,AnltyPbValue,'b*');
hold on;
legend('simulation', 'analytical')
grid on;
title('BER curve vs. SNR (Tilkicioğlu)');
xlabel('SNR(dB)');
ylabel('Pb');
ylim([10^-7, 10^0]);

%----------------------------------------------------------------------

% b) side code   

N = 10^7;                           % 10 million bit size
Eb = 1;                             % average bit energy
Ps1 = 0.25;                         % Probability bit 1
Ps2 = 0.75;                         % Probability bit 0
si = randsrc(1,N,[1,0;Ps1,Ps2]);    % si values
a1 = 1;                             % bit 1 signal value
a2 = -1;                            % bit 0 signal value
ai = 2*si - 1;                      % all ai signal values
max_dB = 17;                        % db value maximum limit

siAllValue = zeros(1,N);            % si vector
simPbValue = zeros(0,max_dB);       % simulation vector
AnltyPbValue = zeros(0,max_dB);     % analytic summaries vector



for db_SNR = 0:max_dB               % loop for all process (0:18 db)
    
   SNR = 10^(db_SNR/10);            % SNR values
   N0 = Eb/SNR;                     % No values
   sigma0 = sqrt(N0);               % sigma value

   z = ai + sigma0*randn(1,N);      % z = (ai + no) equations 

   gamma0 = (N0/2)*log(Ps2/Ps1);    % gamma formula

   siAllValue = z > gamma0;         % decision circuit simulation

   errorVal = si ~= siAllValue;     % Total error for simulation
   Pb = sum(errorVal)/N;            % Pb for simulation
   simPbValue = [simPbValue Pb];    % Pb simulation vector for graph

   % analytic solution 
   analytical_sol = (1-qfunc((gamma0-a1)/sigma0))*Ps1 + (qfunc((gamma0-a2)/sigma0))*Ps2;
   AnltyPbValue = [AnltyPbValue analytical_sol];    % analytic solution 

end

% Solution graphs code for b)

dB_SNR = 0:max_dB;
figure()
semilogy(dB_SNR,simPbValue,'go-');
hold on;
semilogy(dB_SNR,AnltyPbValue,'b*');
hold on;
legend('simulation', 'analytical')
grid on;
title('BER curve vs. SNR (Tilkicioğlu)');
xlabel('SNR(dB)');
ylabel('Pb');
ylim([10^-7, 10^0]);
