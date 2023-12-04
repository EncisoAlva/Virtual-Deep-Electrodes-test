function varargout = process_VirtualChannels( varargin )
% PROCESS_REGION_PRIORS:
% [This function exists for technical purposes]
%
% @========================================================================
% See https://www.overleaf.com/read/cjkvgskgyvjx
% ========================================================================@
%
% Author: Julio Cesar Enciso-Alva, 2023
%         juliocesar.encisoalva@mavs.uta.edu
%
eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
  % Description the process
  sProcess.Comment     = 'Virtual Deep Electrodes';
  sProcess.Category    = 'File';
  sProcess.SubGroup    = 'Sources';
  sProcess.Index       = 1000;
  sProcess.FileTag     = '';
  sProcess.Description = 'https://www.overleaf.com/read/cjkvgskgyvjx';
  % Definition of the input accepted by this process
  sProcess.InputTypes  = {'results','data','raw'};
  sProcess.OutputTypes = {'data','raw'};
  sProcess.nInputs     = 1;
  %sProcess.nMinFiles   = 1;
  %
  sProcess.options.ScoutRadius.Comment = 'Scout Radius: ';
  sProcess.options.ScoutRadius.Type    = 'value';
  sProcess.options.ScoutRadius.Value   = {4, 'mm', 2};   % {Default value, units, precision}
  sProcess.options.ScoutRadius.Class   = 'Debug';
  %
  sProcess.options.magn.Comment   = 'Compute magnitude';
  sProcess.options.magn.Type      = 'checkbox';
  sProcess.options.magn.Value     = 1;
  %
  sProcess.options.debug3.Comment   = 'Debug: Keep 3 coordinates (Unconstrained)';
  sProcess.options.debug3.Type      = 'checkbox';
  sProcess.options.debug3.Value     = 1;
  %
  sProcess.options.debug4.Comment   = 'Process time in batches';
  sProcess.options.debug4.Type      = 'checkbox';
  sProcess.options.debug4.Value     = 1;
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sProcess.Comment;
end

function OutputFiles = Run(sProcess, sInputs)
  % Initialize returned list of files
  OutputFiles = {};

  % ===== GET OPTIONS =====
  % General
  debugFlag   = sProcess.options.debug3.Value;
  magn        = sProcess.options.magn.Value;
  ModalityNew = 'SEEG';
  ModalityDat = 'ECOG';
  Distance    = (sProcess.options.ScoutRadius.Value{1})/1000; % convert to m

  % Get unique channel files 
  AllChannelFiles = unique({sInputs.ChannelFile});
  % ===== LOOP ON FOLDERS =====
  for iChanFile = 1:length(AllChannelFiles)
    % Load channel file
    ChannelMat = in_bst_channel(AllChannelFiles{iChanFile});
    % Get selected sensors
    iChannelsDat = channel_find(ChannelMat.Channel, ModalityDat);
    iChannelsNew = channel_find(ChannelMat.Channel, ModalityNew);
    % Verifications of missing files
    if isempty(iChannelsDat)
      bst_report('Error', sProcess, sInputs, ['Channels "' ModalityDat '" not found in channel file.']);
      return;
    end
    if isempty(iChannelsNew)
      bst_report('Error', sProcess, sInputs, ['Channels "' ModalityNew '" not found in channel file.']);
      return;
    end
    % ===== LOOP ON DATA FILES =====
    % Get data files for this channel file
    iChanInputs = find(ismember({sInputs.ChannelFile}, AllChannelFiles{iChanFile}));
    % Loop on data files
    for iInput = 1:length(iChanInputs)
      % === LOAD DATA ===
      % Load data
      DataFile   = sInputs(iChanInputs(iInput)).FileName;
      ResultsMat = in_bst_results(DataFile);
      iStudyData = sInputs(iChanInputs(iInput)).iStudy;
      % Remove bad channels
      iBadChan   = find(ResultsMat.ChannelFlag == -1);
      iChannelsData = setdiff(iChannelsDat, iBadChan);
      %iChannelsNew
      % Error: All channels tagged as bad
      if isempty(iChannelsData)
        bst_report('Error', sProcess, sInputs, 'All the selected channels are tagged as bad.');
        return;
      end
      % Read locations of intended new channels
      chanPos  = {ChannelMat.Channel.Loc};
      chanType = {ChannelMat.Channel.Type};
      chanName = {ChannelMat.Channel.Name};
      chanUsed = {};
      NewPos = [];
      for i = 1:length(chanName)
        if strcmp(chanType{i}, ModalityNew)
          NewPos = [NewPos; chanPos{i}'];
          chanUsed{end+1} = strcat(chanName{i},'_est');
        end
      end
      % Auto detect if constrined model was used
      if (ResultsMat.nComponents == 1)
        constr = true;
      elseif (ResultsMat.nComponents == 3)
        constr = false;
      else
        bst_report('Error', sProcess, sInputs, 'Cannot run this process on mixed source models.');
        return;
      end
      % If results in compact/kernel mode
      if isfield(ResultsMat, 'ImageGridAmp') && isempty(ResultsMat.ImageGridAmp) && ~isempty(ResultsMat.ImagingKernel)
        % No data file
        if isempty(ResultsMat.DataFile)
          bst_report('Error', sProcess, sInputs, 'No data file defined for this inverse kernel.');
          return;
        end
        % Load data
        DataFileFull = file_fullpath(ResultsMat.DataFile);
        if file_exist(DataFileFull)
          DataMat = load(DataFileFull, 'F', 'Time');
          Time = DataMat.Time;
        else
          bst_report('Error', sProcess, sInputs, 'Data file not found for this inverse kernel.');
          return;
        end
      else
        DataMat = [];
        Time = ResultsMat.Time;
      end
      % === OUTPUT STRUCTURE ===
      % Create structure
      newMat = db_template('matrixmat');
      newMat.Value       = [];
      newMat.ChannelFlag = ones(length(chanUsed),1);
      newMat.Time = Time;
      % If the number of averaged files is defined: use it
      if isfield(ResultsMat, 'nAvg') && ~isempty(ResultsMat.nAvg) % Probably skip
        newMat.nAvg = ResultsMat.nAvg;
      else
        newMat.nAvg = 1;
      end
      if isfield(ResultsMat, 'Leff') && ~isempty(ResultsMat.Leff)
        newMat.Leff = ResultsMat.Leff;
      else
        newMat.Leff = 1;
      end
      % Concatenate new values to existing ones
      newMat.Description = cellfun(@(c) [c ' ' newMat.Comment], chanUsed, 'UniformOutput', false);
      %if ~isempty(scoutStd)
      %  newMat.Std = scoutStd;
      %end
      % Save units, verified consistent if concatenating
      newMat.DisplayUnits = ResultsMat.DisplayUnits;

      % ===== CALCULATE SCOUTS VALUES =====
      [newMat.Value, newMat.smallKernel, newMat.debugChans] = ...
        EstimateChannels( ...
        ResultsMat.ImagingKernel, ...          % W
        DataMat.F(ResultsMat.GoodChannel, :), ... % Y
        Time, ...                              % time
        Distance, ...                          % dist
        NewPos, ...                            % ElecPos
        ResultsMat.GridLoc, ...                % GridPos
        constr, magn, debugFlag, sProcess.options.debug4.Value);

      % === HISTORY ===
      % Re-use the history of the initial file
      newMat.History = ResultsMat.History;
      % History: process name
      newMat = bst_history('add', newMat, 'process', FormatComment(sProcess));
      % History: File name
      newMat = bst_history('add', newMat, 'process', [' - File: ' sInputs(iInput).FileName]);

      % === ADITIONAL DATA ===
      NewMat.GridLoc     = NewPos;
      NewMat.nComponents = 1;

      % === SAVE FILE ===
      % Comment: forced in the options
      if isfield(sProcess.options, 'Comment') && isfield(sProcess.options.Comment, 'Value') && ~isempty(sProcess.options.Comment.Value)
        newMat.Comment = sProcess.options.Comment.Value;Study
      % Comment: Process default (limit size of scout comment)
      %elseif (length(chanUsed) > 20)
      %  newMat.Comment = [ResultsMat.Comment, ' | ' num2str(length(chanUsed)) ' Estimated Deep Electrodes'];
      %elseif ~isempty(chanUsed)
      %  newMat.Comment = [ResultsMat.Comment, ' | Estimated Deep Electrodes (' chanUsed(2:end) ')'];
      else
        newMat.Comment = [ResultsMat.Comment, ' | Estimated Deep Electrodes'];
      end
      % Save new file in database
      % Output study = input study
      ProtocolInfo = bst_get('ProtocolInfo');
      [sStudy, iiStudy] = bst_get('Study', iStudyData );
      % Output filename
      OutFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'matrix_scout');
      % Save on disk
      bst_save(OutFile, newMat, 'v6');
      % Register in database
      db_add_data(iiStudy, OutFile, newMat);
      % Out to list of output files
      OutputFiles{end+1} = OutFile;

      %iAtlas = find(strcmpi({SurfaceMat.Atlas.Name}, sProcess.options.atlas.Value));
    end
  end
end

%%
% USAGE: x = process_notch('Compute', x, sfreq, FreqList)
function [EstChans, smallKernel, debugChans] = EstimateChannels(...
  W, Y, time, dist, ElecPos, GridPos, ...
  constr, magn, debug, debug2)
% Once the Wiener Kernel is computed, it is used to simulate the
% recordings that could be obtained from depth electrodes at given 
% locations. This data is labeled as recordings from virtual electrodes. 
%
% For each virtual electrode, a scout is created containing all the dipoles 
% located within some given distance of the given location. The magnitudes 
% of these dipoles are averaged over each canonical direction; to obtain a
% one-dimensional rcordeing, etiher the magnitude or the first principal 
% value is computed.
%
% SE : # of Surface Electrodes
% TT : # of Timestamps
% GP : # of Grid Points
% DE : # of Deep Electrodes%
%-------------------------------------------------------------------------
% INPUT
%         W  Wiener Kernel, GPxSE
%         Y  Recordings from surface electrodes, SExTT
%      time  Vector of timestamps, 1xTT
%      dist  Max distance dipole-to-electrode in meters, 1x1
%   ElecPos  Positions of intended depth electrodes, DEx3
%   GridPos  Positions of distributed dipoles, GPx3
%
%-------------------------------------------------------------------------
% INPUT (OPTIONAL)
%    constr  Constrained optimization, default = F
%            Fixed orientation (1 column/dipole) or not (3 col/dipole)
%     magn   Magnitude as final result, default = T
%            If  constrained: T=magnitude, F=PCA 1st coordinate
%            If ~constrained: T=abs val, F=no changes
%    debug   If true, and if constrained optimization, give also magnitudes
%
%-------------------------------------------------------------------------
% OUTPUT
%    EstChans  Estimated channels, DExTT
% smallKernel  ...
%
%-------------------------------------------------------------------------
% Author: Julio Cesar Enciso-Alva, 2023
%         juliocesar.encisoalva@uta.edu
%

% === METADATA ===
meta    = [];
meta.SE = size(Y, 1);      % Surface Electrodes
meta.TT = size(Y, 2);      % Timestamps
meta.GP = size(GridPos,1); % Grid Points
meta.DE = size(ElecPos,1); % Deep Electrodes

% === SIZE VERIFICATION ===
if meta.SE ~= size(W, 2)
  disp('ERROR: Size mismatch: number of Deep Electrodes.')
  return
end
if constr
  if meta.GP ~= size(W, 1)
    disp('ERROR: Size mismatch: number of Grid Points.')
     return
  end
else
  if meta.GP*3 ~= size(W, 1)
    disp('ERROR: Size mismatch: number of Grid Points.')
    return
  end
end
if size(time,2) ~= meta.TT
  disp('ERROR: Size mismatch: number of Timestamps.')
  return
end
if size(GridPos,2)~=3 || size(ElecPos,2)~=3
  disp('ERROR: Size mismatch: space coordinates must be 3.')
  return
end

% === PROCESSING ===
% create scouts, nearby dipoles
scout = cell( meta.DE, 1 );
idx   = 1:meta.GP;
for ii = 1:meta.DE
  scout{ii} = idx(vecnorm( GridPos-ElecPos(ii,:), 2, 2 ) < dist);
end

% average kernel within scouts
if constr
  shortW = zeros(  meta.DE, meta.SE);
else
  shortW = zeros(3*meta.DE, meta.SE);
end
for ii = 1:meta.DE
  if constr
    shortW(ii,:) = mean( W(scout{ii},:), 1 );
  else
    disp(ii)
    shortW(3*(ii-1)+1,:) = mean( W(3*(scout{ii}-1)+1,:), 1 );
    shortW(3*(ii-1)+2,:) = mean( W(3*(scout{ii}-1)+2,:), 1 );
    shortW(3*(ii-1)+3,:) = mean( W(3*(scout{ii}-1)+3,:), 1 );
  end
end
if debug
  smallKernel = shortW;
else
  smallKernel = [];
end

% construction of estimated channels
rawChans = shortW*Y;
if debug
  debugChans = rawChans;
else
  debugChans = [];
end
if constr
  if magn
    EstChans = abs(rawChans);
  else
    EstChans = rawChans;
  end
else
  % making it a tensor for ease of notation
  idx = 0:(meta.DE-1);
  tensChans = zeros(meta.DE,meta.TT,3);
  tensChans(:,:,1) = rawChans(3*idx+1,:);
  tensChans(:,:,2) = rawChans(3*idx+2,:);
  tensChans(:,:,3) = rawChans(3*idx+3,:);
  if magn
    EstChans = vecnorm(tensChans,2,3);
  else
    EstChans = zeros(meta.DE,meta.TT);
    if debug2
      n_batches  = ceil(meta.TT/5000);
      batch_size = floor(meta.TT/n_batches);
    else
      n_batches  = 1;
      batch_size = meta.TT;
    end
    for ii = 1:meta.DE
      for batch = 1:n_batches
        if batch < n_batches
          b_time = ((batch-1)*batch_size+1):(batch*batch_size);
        else
          b_time = ((batch-1)*batch_size+1):meta.TT;
        end
        ChanVec = reshape( tensChans(ii,b_time,:), max(size(b_time)), 3);
        [EIG,~] = eigs(cov( ChanVec' ));
        EstChans(ii,b_time) = EIG(:, 1)';
      end
    end
  end
end

end