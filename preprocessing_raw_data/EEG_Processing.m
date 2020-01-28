% preprocessing and segmentation of the raw EEG data
% EEGLab Toolbox is needed 

% First, the raw data is processed with a Butterworth fourth order bandpass filter of 1-30 Hz. Artifacts above and below 100 microvolts are removed.
% Each electrode is re-referenced by the average reference.
% Secondly, the filtered EEG data is segmented as follows.
%	  (i) The data is divided in two parts, each part corresponding to the blocks QIT and TIQ.
%	 (ii) From each block, two segments are extracted, one for each chain, Quaternary and Ternary. 
%	(iii) The EEG chunks corresponding to the Ternary/Quaternary chain are segmented as follows. For each symbol {0,1,2}, EEG segments of 450ms
%         are extracted from 50ms before the stimulus input to 400ms after it. Then, a baseline correction is performed using the signal collected
%         50 ms before each event input. 


%Author : Aline Duarte, Noslen Hernandez
%Date   : 01/2020

clear all

% directory to load the raw data
base_rootdir = '/Users/numec/Desktop/Raw_Data';
% directory to save the preprocessed data
base_destdir = '/Users/numec/Google Drive/Mija_files/Final_preprocessing';

% parameters used in the segmentation and base line correction
tx_amost = 250;
t_min = -50;
t_max = 400;
rmbase_min = -50; 
rmbase_max = 0;

% if the directory to save the preprocessed does not exist, create it
if ~exist( base_destdir, 'dir' )
    mkdir( base_destdir );
end
rootdir = fullfile( base_rootdir);
files = dir( fullfile( rootdir, '*.raw' ) );

% for each volunteer
for v = [ 1:3  5:length(files) ]	% we skip volunteer V04 as indicate in the article 
    
	EEG = pop_readegi( fullfile( rootdir, files(v).name ), [], [], 'GSN_HydroCel_129.sfp' );
    EEG.setname = files(v).name(1:end-4);
    EEG.filename = EEG.setname;
    volname = EEG.filename(1:3);
    
	% filtering
    HP = 1;
    LP = 30;
    ddata = double(EEG.data);
    order = 4;
    datafiltered = filterpass( ddata', HP,  EEG.srate, 'high',  order);
    datafiltered = filterpass( datafiltered, LP,  EEG.srate, 'low' , order);
    EEG.data = datafiltered';
    
	% channel rejection and re-reference
    std_fac = [-3 3];
    freqs = [1:29 59; 2:30 61];
    freqs = reshape(freqs,2,size(freqs,2));
    [EEG_chan_rej indelec] = pop_rejchanspec( EEG, 'freqlims', freqs, 'stdthresh', repmat( std_fac, size(freqs,1), 1 ) );
    EEG_chan_rej = pop_reref( EEG_chan_rej,[] );
    EEG = EEG_chan_rej;
    
	% take only the 10_20 channel set
    set_10_20 = {'E9','E11','E22','E24','E33','E36','E45','E52','E58','E62','E70','E83','E92', 'E96','E104','E108','E122','E124'};         
    EEG = pop_select( EEG, 'channel', set_10_20 ); 
	
	% sharing in blocks
    Event = squeeze( struct2cell( EEG.event ) );
    foundevents = find( strcmp( Event(1,:), 'sil ' ) ); % vector of positions in which the event appears in eventcells
    idx = [];
    for s = 1 : (length(foundevents)-1)
        if foundevents(s+1)-foundevents(s) == 134 % number of stimuli in one minute of rhythms 
            idx = [idx, foundevents(s)];
        end
    end
    indices = cell2mat( {EEG.event(idx).latency} ); % vector of latencies of that event
    
	%
     if length( indices ) == 18
        EEG_b1 = pop_select( EEG, 'point', [indices(1)  cell2mat({EEG.event(idx(9)+134).latency}) + t_max] );
        EEG_b1.filename = [ EEG_b1.filename '_B1'];
        EEG_b2 = pop_select( EEG, 'point', [indices(10)  cell2mat({EEG.event(idx(18)+134).latency}) + t_max] );
        EEG_b2.filename = [EEG_b2.filename '_B2' ];
    else
        warning( sprintf( 'could not segment file %s, it is not an 18min sample', files(v).name ) )
    end
    
	% sharing rhythms
    inds = strfind( EEG.setname, '_' );
    ALLEEG_aux = {};
    for r = 1:3
        rhythm = EEG.setname(inds + r);
        if strcmpi( rhythm , 'T' ),
           rhythm_name = 'Ter'; ind = 1;
           elseif strcmpi( rhythm, 'Q' );
           rhythm_name = 'Qua'; ind = 2;
           elseif strcmpi( rhythm, 'I' );
           rhythm_name = 'Ind'; ind = 3;
        end
        EEG_rhythm = pop_select( EEG, 'point', [indices(3*(r-1)+1)  cell2mat( {EEG.event(idx(3*r)+134).latency}) + t_max ] );
        EEG_rhythm.filename = [EEG_rhythm.setname(1:end-4) '_' rhythm_name '-b1' ];
        ALLEEG_aux = eeg_store(ALLEEG_aux, EEG_rhythm, ind);
        ALLEEG_aux(ind).setname = [ rhythm_name '_'  EEG.setname];
    end

    for r = 4:6
        rhythm = EEG.setname(end-r+4);
        if strcmpi( rhythm , 'T' ),
          rhythm_name = 'Ter'; ind = 1;
        elseif strcmpi( rhythm, 'Q' );
          rhythm_name = 'Qua'; ind = 2;
        elseif strcmpi( rhythm, 'I' );
          rhythm_name = 'Ind'; ind = 3;
        end
        EEG_rhythm = pop_select( EEG, 'point', [indices(3*(r-1)+1)  cell2mat( {EEG.event(idx(3*r)+134).latency}) + t_max ] );
        EEG_rhythm.filename = [EEG_rhythm.setname(1:end-4) '_' rhythm_name '-b2' ];
        ALLEEG_aux = eeg_store(ALLEEG_aux, EEG_rhythm, ind+3);
        ALLEEG_aux(ind+3).setname = [ rhythm_name '_' EEG.setname];
    end
	
    % merging dataset
    Event_ter = squeeze( struct2cell( ALLEEG_aux(1).event ) );
    foundev = find( strcmp( Event_ter(1,:), 'sil ' )); % the same indexes for the 6 files
 
    ALLEEG_aux(1).event( [foundev, foundev(2:end-1) + 1] ) = [];
    ALLEEG_aux(4).event( [foundev, foundev(2:end-1) + 1, 2] ) = [];
    EEG_ter = pop_mergeset(ALLEEG_aux(1), ALLEEG_aux(4));
    Event_ter = squeeze( struct2cell( EEG_ter.event ) );
    erase = find( strcmp( Event_ter(1,:), 'boundary' ));
    EEG_ter.event(erase) = [];
    EEG_ter.setname = ALLEEG_aux(1,1).setname;
    
    ALLEEG_aux(2).event( [foundev, foundev(2:end-1) + 1] ) = [];
    ALLEEG_aux(5).event( [foundev, foundev(2:end-1) + 1,2] ) = [];
    EEG_qua = pop_mergeset(ALLEEG_aux(2), ALLEEG_aux(5));
    Event_qua = squeeze( struct2cell( EEG_qua.event ) );
    erase = find( strcmp( Event_qua(1,:), 'boundary' ));
    EEG_qua.event(erase) = [];
    EEG_qua.setname = ALLEEG_aux(1,1).setname;
    
    ALLEEG_aux(3).event([foundev, foundev(2:end-1)+1]) = [];
    ALLEEG_aux(6).event([foundev, foundev(2:end-1)+1,2]) = [];
    EEG_ind = pop_mergeset(ALLEEG_aux(3),ALLEEG_aux(6));
    Event_ind = squeeze( struct2cell( EEG_ind.event ) );
    erase = find( strcmp( Event_ind(1,:), 'boundary' ));
    EEG_ind.event(erase) = [];
    EEG_ind.setname = ALLEEG_aux(1,1).setname;
    
    %
    Event = squeeze( struct2cell( EEG_ter.event ) );
    X_ter(2,:) = cell2mat(Event(2,:));
    for i = 1 : length(Event)
        if  strcmp( Event(1,i), 'v2  ' )
            X_ter(1,i) = str2num('2');
        elseif strcmp( Event(1,i), 'v1a ' )
            X_ter(1,i) = 1;
        elseif strcmp( Event(1,i), 'v1b ' )
            X_ter(1,i) = 1;
        elseif strcmp( Event(1,i), 'miss' )
            X_ter(1,i) = 0;
        end
    end
    data.X_ter = X_ter(1,:);
    
    Event = squeeze( struct2cell( EEG_qua.event ) );
    X_qua(2,:) = cell2mat(Event(2,:));
    for i = 1 : length(Event)
        if  strcmp( Event(1,i), 'v2  ' )
            X_qua(1,i) = str2num('2');
        elseif strcmp( Event(1,i), 'v1a ' )
            X_qua(1,i) = 1;
        elseif strcmp( Event(1,i), 'v1b ' )
            X_qua(1,i) = 1;
        elseif strcmp( Event(1,i), 'v0  ' )
            X_qua(1,i) = 0;
        elseif strcmp( Event(1,i), 'miss' )
            X_qua(1,i) = 0;
        end
    end
    data.X_qua = X_qua(1,:);
    
    Y_ter = {};
    Y_qua = {};
    for e = 1 : EEG.nbchan
        Y_ter{1,e} = EEG.chanlocs(e).labels;
        Y_ter_e = EEG_ter.data(e,:);
        Y_ter{2,e} = segmentation( X_ter, Y_ter_e, tx_amost, t_min, t_max, rmbase_min, rmbase_max  );
        
        Y_qua{1,e} = EEG.chanlocs(e).labels;
        Y_qua_e = EEG_qua.data(e,:);
        Y_qua{2,e} = segmentation( X_qua, Y_qua_e, tx_amost, t_min, t_max, rmbase_min, rmbase_max  );
    end
    data.Y_ter=Y_ter;
    data.Y_qua=Y_qua;
    
	% save the preprocessed data
    filename = fullfile(base_destdir , EEG.setname(1:3)); 
    save(filename, 'data' )
end
