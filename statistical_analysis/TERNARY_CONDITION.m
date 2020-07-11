% Statistical Analysis for the Ternary condition

% number of the electrodes used in the statistical analysis
target_electrodes = [9, 11, 22, 24, 33, 36, 45, 52, 58, 62, 70, 83, 92, 96, 104, 108, 122, 124];

% alphabet
A = [0,1,2];

% values of the parameters
height_complete_tree = 3;	% height of the complete tree
n_BM = 5000;				% number of Brownian motions
alpha = 0.05;				% significance of the ks-test
beta = 0.05;				% Bernoulli parameter

% variable to store the estimated context tree for each participant and each electrode
tree_array_ter = cell(19,18);

% for each participant and each electrode, estimate a context tree from the EEG data
participant_counter = 1;

for s = [1:3 5:20] 		% for each participant

    % load the data of the subject
    name = ['V' num2str(s,'%02d') '.mat'];
    load(['preprocessed_data/' name]);
    
	% load the sequence of stimuli
    X = data.X_ter;
    
    % verify if electrode information match
    electrodes = cellfun(@(x) str2double(x(2:end)), data.Y_ter(1,:), 'UniformOutput', true); 
    
    if sum(electrodes == target_electrodes) == 18
        for electrode_counter = 1 : 18
            disp(['Processing Participant ' name(1:end-4) ', electrode ' num2str(electrode_counter) '...']);
            % load the EEG data
			Y = data.Y_ter{2,electrode_counter};     
			% fix the seed to use the same Brownian motions for all participants and electrodes
            rng(1)
			% estimate the context tree	
            tree = estimate_functionalSeqRoCTM(X, Y, A, height_complete_tree, n_BM, alpha, beta);
            % store the context tree
			tree_array_ter{participant_counter, electrode_counter} = tree; 
        end
    else
        disp('EEG information about electrodes does not match the selected electrodes.')
    end
    participant_counter = participant_counter + 1;
end

save('ternary_trees','tree_array_ter');

% compute the mode context tree on each electrode
mode_context_tree_ter = cell(18,1);     % mode context tree for each electrode

for e = 1 : 18
    trees = tree_array_ter(:,e);
	mode_context_tree_ter{e} = mode_tree(trees, A);
end
