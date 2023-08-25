% Inputs modulation: Vector of source and interference modulation formats.
%                    Accepted values: "QAM", "PSK", and "WGN".
%        Mmod:       Vector of modulation orders of the modulation vector.
%                    Note that, for "WGN", the modulation order is
%                    unnecessary.
%        f:          Operation frequency in Hz.
%        phi:        Angle of the transmmiters (azimuth) in degrees.
%        theta:      Angle of the transmitters (zenith) in degrees.
%        desired:    Vector to select sources and interferences. If 1, the
%                    respective signal is from a source, if 0 it is from an
%                    interference.
%        lenData:    Number of samples to the dataset.
%        SINRdB:     Signal-to-interference plus noise ratio in dB.
%        SNRdBs:     Signal-to-noise ratio of sources in dB.
%        SNRdBi:     Signal-to-noise ratio of interferences in dB.

% Output SetIn:      Input dataset. It is a matrix with 6 rows (due to the 
%                    six dipoles at the Rx) and lenData columns.
%        SetOut:     Output dataset. It is a matrix with a number of rows
%                    given by the number of users and lenData columns.

% Example of input parameters:
% f          = 850e6;                                                                                          
% SINRdB     = 20;                                                                 
% SNRdBs     = 25;                                                                   
% SNRdBi     = 20;                                                                                                                                      
% phi        = [ 1  60  90  120  160  200  240  260  280  300  330];
% theta      = [90  90  90   90   90   90   90   90   90   90   90];
% desired    = [ 1   0   0    1    0    0    1    0    0    0    0];
% modulation = ["QAM" "WGN" "QAM" "PSK" "QAM" "WGN" "QAM" "WGN" "QAM" ...
%               "PSK" "PSK"];
% Mmod       = [ 4  0  64  8  256  0  16  0  64  16  8];
% lenData    = 3e4;

% Created by Kayol Soares Mayer
% Last updated 2023/08/25

% Cite:
% @ARTICLE{Mayer2022,
%  author={Mayer, Kayol Soares and M\"uller, Candice and Soares, 
%          Jonathan Aguiar and de Castro, Fernando C\'esar Comparsi and 
%           Arantes, Dalton Soares},
%  journal={IEEE Wireless Communications Letters}, 
%  title={{Deep phase-transmittance RBF neural network for beamforming with 
%          multiple users}}, 
%  year={2022},
%  volume = {11},
%  number = {7},
%  pages = {1498--1502},
%  }

function [SetIn,SetOut] = create_data_set_beam(modulation,Mmod,f,phi,...
                                theta,desired,lenData,SINRdB,SNRdBs,SNRdBi)

    % Linear SNR of sources
    SNRs = 10^(SNRdBs/10);
    % Linear SNR of interferences
    SNRi = 10^(SNRdBi/10);
    % Linear SINR
    SINR = 10^(SINRdB/10);
    
    % Number of sources
    Ns = length(find(desired));
    % Number of interferences
    Ni = length(find(~desired));
    
    % Normalization factor to the desired SINR
    iota = (Ni/Ns)*(1/SNRi+1)/(1/SINR-1/SNRs);

    % Standard deviations for AWGN noises of sources and interferences
    StdDevS  = sqrt((1/SNRs)/2);
    StdDevI  = sqrt((1/(iota*SNRi))/2);
    
    % Vector of AWGN standard deviations
    StdDevVect = (desired*StdDevS+~desired*StdDevI)';
    
    % speed of light in vacuum [m/s]
    c=299792458;
    % wavelength
    lambda=c/f;
    % propagation constant
    beta=2*pi/lambda;
    % dipoles length
    L=0.5*lambda;
    % dipoles spacing
    s=0.25*lambda;
    % dipoles cooordinates
    coord=[     s          0         0;
            s*cosd(60) s*sind(60)    0;
           -s*cosd(60) s*sind(60)    0;
               -s          0         0;
           -s*cosd(60) -s*sind(60)   0;
            s*cosd(60) -s*sind(60)   0]';  
        
    % Matrix of self and mutual impedances
    Z = [78.424+45.545i 46.9401-32.6392i -0.791-41.3825i -14.4422-34.4374i -0.791-41.3825i 46.9401-32.6392i;...    
         46.9401-32.6392i 78.424+45.545i 46.9401-32.6392i -0.791-41.3825i -14.4422-34.4374i -0.791-41.3825i;...
         -0.791-41.3825i 46.9401-32.6392i 78.424+45.545i 46.9401-32.6392i -0.791-41.3825i -14.4422-34.4374i;...
         -14.4422-34.4374i -0.791-41.3825i 46.9401-32.6392i 78.424+45.545i 46.9401-32.6392i -0.791-41.3825i;...
         -0.791-41.3825i -14.4422-34.4374i -0.791-41.3825i 46.9401-32.6392i 78.424+45.545i 46.9401-32.6392i;...
         46.9401-32.6392i -0.791-41.3825i -14.4422-34.4374i -0.791-41.3825i 46.9401-32.6392i 78.424+45.545i];

    % Load impedance
    ZT = 50; %[ohms]

    % Dipoles self impedance
    ZA = Z(1,1);

    % Coupling matrix
    C = (ZT+ZA)*inv((Z+ZT*eye(6)));
    
    % Matrix of relative intensity of Etheta
    Xm=diag((lambda/(pi*sin(beta*L/2)))*((cos(L*pi*cosd(theta')/lambda)-cos(pi*L/lambda))./sind(theta')));
    
    % Matrix of Tx angular positions
    Omega = [sind(theta)'.*cosd(phi)' sind(theta)'.*sind(phi)' cosd(theta)'];
    
    Pi = exp(-2i*pi*(Omega*coord)/lambda);
    
    % Steering vectors
    psi = Pi.'*Xm;
    
    % Create symbols of sources and interferences
    SetOut = zeros(Ns+Ni,lenData);
    for ii=1:1:Ns+Ni
        % QAM symbols
        if modulation(ii) == "QAM"
            if Mmod(ii) ~= 0
                SetOut(ii,:) =   (((randi(sqrt(Mmod(ii)),1,lenData)-1)*2)-(sqrt(Mmod(ii))-1))+...
                              1i*(((randi(sqrt(Mmod(ii)),1,lenData)-1)*2)-(sqrt(Mmod(ii))-1));
            else
                SetOut(ii,:) =   randn(1,lenData)+1i*randn(1,lenData);
            end
        % PSK symbols
        elseif modulation(ii) == "PSK"     
            if Mmod(ii) ~= 0
                pskAng = randi(Mmod(ii),1,lenData)*2*pi/Mmod(ii);
                SetOut(ii,:) =   cos(pskAng)+1i*sin(pskAng);
            else
                SetOut(ii,:) =   randn(1,lenData)+1i*randn(1,lenData);
            end
        % WGN noise
        else
            SetOut(ii,:) =   randn(1,lenData)+1i*randn(1,lenData);
        end
    end
    
    % Compute the source and interference powers to normalization
    P = sum((abs(SetOut-mean(SetOut,2))).^2,2)/size(SetOut,2);
    
    % Normalize the powers to 1
    SetOut = SetOut./sqrt(P);
   
    % Interference normalizations to the desired SINR
    SetOut(find(~desired),:) = SetOut(find(~desired),:)/sqrt(iota);
    
    % Create the data that impinged on the beamforming
    SetIn = (C*psi)*(SetOut+StdDevVect.*(randn(size(SetOut))+1i*randn(size(SetOut))));
    
    % Compute the signal powers in each RX dipole
    P=sum((abs(SetIn-mean(SetIn,2))).^2,2)/size(SetIn,2);
    
    % Normalize the signals in each dipole to have a unitary power
    SetIn = SetIn./sqrt(P);
    
    % Add nonlinearities
    SetIn = SetIn - 0.1*SetIn.^3 - 0.05*SetIn.^5;
    
    % Select desired outputs
    SetOut = SetOut(find(desired),:);
    
end

