function S = segmentation(X, Y, tx_amost, t_min, t_max, rmbase_min, rmbase_max)
% SEGMENTATION segment the EEG signal and do baseline correction
%
% INPUTS
%
% X 		 : sequence of stimuli containing in the second row the time at which each stimuli was presented
% Y 		 : EEG signal
% tx_amost	 : sample rate
% t_min 	 : initial time for segmentation
% t_max 	 : final time for segmentation
% rmbase_min : initial time for baseline correction
% rmbase_max : final time for baseline correction
%
% OUTPUTS
%
% S 		 : EEG segments

%Author : Aline Duarte, Noslen Hernandez
%Date   : 01/2020

nHz = 1000; %unit of frequency 

nb_inf = ceil( t_min*tx_amost/nHz );
nb_sup = ceil( t_max*tx_amost/nHz );
rm_inf = ceil( rmbase_min*tx_amost/nHz );
rm_sup = ceil( rmbase_max*tx_amost/nHz );

S = zeros(abs(nb_inf) + abs(nb_sup+1), length(X));

%
for k = 1 : length(X) 
    I = Y(X(2,k) + nb_inf : X(2,k) + nb_sup);
    rm_mean = mean(Y(X(2,k) + rm_inf : X(2,k) + rm_sup));

    S(:,k) = I - rm_mean;
end

end

