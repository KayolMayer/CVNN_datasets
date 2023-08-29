% Code to create an equalization dataset for CVNNs

% Inputs lenData: Number of samples to the dataset.

% Output SetIn:   Input dataset. It is a vector with lenData columns.
%        SetOut:  Output dataset. It is a vector with lenData columns.

% Example of input parameters:
% lenData = 5e3;

% Created by Kayol Soares Mayer
% Last updated 2023/08/29

% Cite:
%@ARTICLE{Cha1995,
%  author={Inhyok Cha and Kassam, S.A.},
%  journal={IEEE Journal on Selected Areas in Communications}, 
%  title={{Channel equalization using adaptive complex radial basis function networks}}, 
%  year={1995},
%  volume={13},
%  number={1},
%  pages={122--131},
%}

function [SetIn,SetOut] = create_dataset_equalization(lenData)
    
    % Initialize the input dataset with zeros to speed up the code.
    SetIn  = zeros(1,lenData);
    % Create the output dataset with QPSK symbols
    SetOut = (randi([0,1],1,lenData)*2-1)/sqrt(2)+...
                                   1i*(randi([0,1],1,lenData)*2-1)/sqrt(2);
    
    % Initialize the channel delays
    a1 = -1/sqrt(2)-1i/sqrt(2);
    a2 = -1/sqrt(2)+1i/sqrt(2);

    % Loop through the dataset samples
    for ii = 1:1:lenData
        % Get the current data sample from QPSK symbols
        a0 = SetOut(ii);
        
        % Apply the channel convolution
        aux = (0.34-1i*0.27)*a0+(0.87+1i*0.43)*a1+(0.34-1i*0.21)*a2;
        
        % Apply nonlinearities and add AWGN noise
        SetIn(ii) = aux + 0.1*aux^2+ 0.05*aux^3+...
                           sqrt(0.01)*(randn()/sqrt(2)+1i*randn()/sqrt(2));
        
        % shift the channel samples for the next channel operations
        a2=a1;
        a1=a0;

    end
end

