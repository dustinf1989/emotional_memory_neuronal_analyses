function [spike] = ft_read_spike(filename, varargin) % Edited to remove unit and waveform fields

% FT_READ_SPIKE reads spike timestamps and waveforms from various data formats.
%
% Use as
%  [spike] = ft_read_spike(filename, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'spikeformat' = string, described the fileformat (default is automatic)
%
% The following file formats are supported
%   'blackrock_nev'
%   'mclust_t'
%   'neuralynx_ncs'
%   'neuralynx_nse'
%   'neuralynx_nst'
%   'neuralynx_ntt'
%   'neuralynx_nts'
%   'neuroshare'
%   'neurosim_spikes'
%   'nwb'
%   'plexon_ddt'
%   'plexon_nex'
%   'plexon_nex5'
%   'plexon_plx'
%   'wave_clus'
%
% The output spike structure usually contains
%   spike.label     = 1xNchans cell-array, with channel labels
%   spike.waveform  = 1xNchans cell-array, each element contains a matrix (Nleads x Nsamples X Nspikes)
%   spike.waveformdimord = '{chan}_lead_time_spike'
%   spike.timestamp = 1xNchans cell-array, each element contains a vector (1 X Nspikes)
%   spike.unit      = 1xNchans cell-array, each element contains a vector (1 X Nspikes)
% and is described in more detail in FT_DATATYPE_SPIKE
%
% See also FT_DATATYPE_SPIKE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT

% Copyright (C) 2007-2021 Robert Oostenveld, Arjen Stolk
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% optionally get the data from the URL and make a temporary local copy
filename = fetch_url(filename);

if ~exist(filename,'file')
  ft_error('File or directory does not exist')
end

% get the options
spikeformat = ft_getopt(varargin, 'spikeformat');

% optionally get the data from the URL and make a temporary local copy
filename = fetch_url(filename);

if isempty(spikeformat)
  % only do the autodetection if the format was not specified
  spikeformat = ft_filetype(filename);
end

switch spikeformat
  case {'neurosim_spikes' 'neurosim_ds'}
    spike = read_neurosim_spikes(filename);
    
  case {'neuralynx_ncs' 'plexon_ddt'}
    % these files only contain continuous data
    ft_error('file does not contain spike timestamps or waveforms');
    
  case 'matlab'
    % plain MATLAB file with a single variable in it
    load(filename, 'spike');
    
  case 'mclust_t'
    fp = fopen(filename, 'rb', 'ieee-le');
    H = ReadHeader(fp);
    fclose(fp);
    % read only from one file
    S = read_mclust_t({filename});
    spike.hdr = H(:);
    [p, f, x] = fileparts(filename);
    spike.label     = {f};  % use the filename as label for the spike channel
    spike.timestamp = S;
%     spike.waveform  = {};   % this is unknown
%     spike.unit      = {};   % this is unknown
    spike.hdr       = H;
    
  case 'wave_clus'
    load(filename, 'cluster_class', 'spikes', 'par'); % load the mat file
    clusters = sort(unique(cluster_class(:,1))); % detected clusters
    clusters(clusters==0) = []; % remove rejected cluster (indexed by zeros)
    nclust = numel(clusters);
    [p, f, x] = fileparts(filename);
    t = tokenize(f, '_'); % extract channel name
    spike.label     = cell(1,nclust);
    spike.unit      = cell(1,nclust);
    spike.waveform  = cell(1,nclust);
    spike.timestamp = cell(1,nclust);
    spike.hdr       = par;
    for cl = 1:nclust
      unit_idx                  = cluster_class(:,1)==cl;
      spike.label{cl}           = [t{2} '-' num2str(cl)];
      spike.timestamp{cl}       = cluster_class(unit_idx,2)';
      spike.waveform{cl}(1,:,:) = spikes(unit_idx,:)';
      spike.unit{cl}            = cluster_class(unit_idx,1)';
    end
    fprintf('note that wave_clus timestamps are typically expressed in millisec and not in samples\n')

  case 'combinato'
    load(filename, 'cluster_class', 'spikes'); % load the mat file
    clusters = sort(unique(cluster_class(:,1))); % detected clusters
    clusters(clusters==0) = []; % remove rejected cluster (indexed by zeros)
    nclust = numel(clusters);
    [p, f, x] = fileparts(filename);
    t = tokenize(f, '_'); % extract channel name
    spike.label     = cell(1,nclust);
    spike.unit      = cell(1,nclust);
    spike.waveform  = cell(1,nclust);
    spike.timestamp = cell(1,nclust);
%     spike.hdr       = par; % Needed to remove this part
    for cl = 1:nclust
      unit_idx                  = cluster_class(:,1)==cl;
      spike.label{cl}           = [t{2} '-' num2str(cl)];
      spike.timestamp{cl}       = cluster_class(unit_idx,2)';
      spike.waveform{cl}(1,:,:) = spikes(unit_idx,:)';
      spike.unit{cl}            = cluster_class(unit_idx,1)';
    end
    fprintf('note that wave_clus timestamps are typically expressed in millisec and not in samples\n')

   case 'combinato_df'
%     load(filename, 'class_times', 'spikes'); % load the mat file
% spike_class is each small class inside a group. What's the point of keeping them?
    load(filename); % load the mat file
    nonartifacts = group_times(:,2) > 0;
    group_times = group_times(nonartifacts,:);
    spikes = spikes(nonartifacts,:);
    clusters = sort(unique(group_times(:,2))); % detected clusters
%     clusters(spike_group==0) = []; % remove rejected cluster (indexed by zeros)
%     clusters(spike_group==-1) = []; % remove rejected cluster (indexed by zeros)
    nclust = numel(clusters);
    [p, f, x] = fileparts(filename);
    t = tokenize(f, '_'); % extract channel name
    spike.label     = cell(1,nclust);
    spike.unit      = cell(1,nclust);
    spike.waveform  = cell(1,nclust);
    spike.timestamp = cell(1,nclust);
%     spike.hdr       = par; % Needed to remove this part
    for cn = 1:nclust
      cl = clusters(cn);
      unit_idx                  = group_times(:,2)==cl;
      spike.label{cn}           = [t{2} '-' num2str(cl)];
      spike.timestamp{cn}       = group_times(unit_idx,1)';
      spike.waveform{cn}(1,:,:) = spikes(unit_idx,:)';
      spike.unit{cn}            = group_times(unit_idx,2)';
    end
    fprintf('Timestamps from NeuraLynx data are in microseconds\n')

 case 'combinato_pn'
%     load(filename, 'class_times', 'spikes'); % load the mat file
% spike_class is each small class inside a group. What's the point of keeping them?
    load(filename); % load the mat file
    % positive spikes
    nonartifacts = group_times(:,2) > 0;
    group_times = group_times(nonartifacts,:);
    spikes = spikes(nonartifacts,:);
    clusters = sort(unique(group_times(:,2))); % detected clusters
    % negative spikes
    neg_nonartifacts = neg_group_times(:,2) > 0;
    neg_group_times = neg_group_times(neg_nonartifacts,:);
    neg_spikes = neg_spikes(neg_nonartifacts,:);
    neg_clusters = sort(unique(neg_group_times(:,2))); % detected clusters

%     clusters(spike_group==0) = []; % remove rejected cluster (indexed by zeros)
%     clusters(spike_group==-1) = []; % remove rejected cluster (indexed by zeros)
    nclust = numel(clusters);
    nclust_neg = numel(neg_clusters);
    [p, f, x] = fileparts(filename);
    t = tokenize(f, '_'); % extract channel name
    spike.label     = cell(1,nclust+nclust_neg);
    spike.unit      = cell(1,nclust+nclust_neg);
    spike.waveform  = cell(1,nclust+nclust_neg);
    spike.timestamp = cell(1,nclust+nclust_neg);
%     spike.hdr       = par; % Needed to remove this part
    for cn = 1:nclust
      cl = clusters(cn);
      unit_idx                  = group_times(:,2)==cl;
      spike.label{cn}           = [t{4} '_pos_' num2str(cl)];
      spike.timestamp{cn}       = group_times(unit_idx,1)';
      spike.waveform{cn}(1,:,:) = spikes(unit_idx,:)';
      spike.unit{cn}            = group_times(unit_idx,2)';
    end

    for cn = nclust+1:nclust_neg+nclust
      cl = neg_clusters(cn-nclust);
      unit_idx                  = neg_group_times(:,2)==cl;
      spike.label{cn}           = [t{4} '_neg_' num2str(cl)];
      spike.timestamp{cn}       = neg_group_times(unit_idx,1)';
      spike.waveform{cn}(1,:,:) = neg_spikes(unit_idx,:)';
      spike.unit{cn}            = neg_group_times(unit_idx,2)';
    end

    fprintf('Timestamps from NeuraLynx data are in microseconds\n')

   case 'osort'
    load(filename)
    if ~isempty(assignedNegative)
        clusters = sort(unique(assignedNegative)); % detected clusters
        clusters(clusters==0) = []; % remove rejected cluster (indexed by zeros)
        nclust = numel(clusters);
        [p, f, x] = fileparts(filename);
        t1 = tokenize(f, '_'); % extract channel name
        t = [t1{1}, '_', t1{2}];
        % Use regular expression to find the number '1' in the string
        spike.label     = cell(1,nclust);
        spike.unit      = cell(1,nclust);
        spike.waveform  = cell(1,nclust);
        spike.timestamp = cell(1,nclust);
        for cl = 1:nclust
          unit_idx                  = assignedNegative==clusters(cl);
          % newTimestampsNegative(assignedNegative==unit_idx)
          spike.label{cl}           = strcat(t,'_',num2str(clusters(cl))); %clusters(cl)
          spike.timestamp{cl}       = newTimestampsNegative(unit_idx); %class_times(unit_idx,1)';
          spike.waveform{cl}(1,:,:) = newSpikesNegative(unit_idx,:)';
          spike.unit{cl}            = ones(1,length(unit_idx))*cl;
        end
    else
        warning('No Neurons were detected - Maybe there are positive ones?')
    end

   case 'combinato_group'
    load(filename, 'class_times', 'spikes','spike_group'); % load the mat file
    clusters = sort(unique(spike_group)); % detected clusters
    clusters(clusters==0) = []; % remove rejected cluster (indexed by zeros)
    nclust = numel(clusters);
    [p, f, x] = fileparts(filename);
    t = tokenize(f, '_'); % extract channel name
    spike.label     = cell(1,nclust);
    spike.unit      = cell(1,nclust);
    spike.waveform  = cell(1,nclust);
    spike.timestamp = cell(1,nclust);
%     spike.hdr       = par; % Needed to remove this part
    for cl = 1:nclust
      unit_idx                  = spike_group==cl;
      spike.label{cl}           = [t{2} '-' num2str(cl)];
      spike.timestamp{cl}       = class_times(unit_idx,1)';
      spike.waveform{cl}(1,:,:) = spikes(unit_idx,:)';
      spike.unit{cl}            = spike_group';
    end
    fprintf('note that combinato timestamps are typically expressed in microsec and not in samples\n')

  case 'neuralynx_ntt'
    % single channel stereotrode file, read all records
    ntt = read_neuralynx_ntt(filename);
    clusters = sort(unique(ntt.CellNumber)); % detected clusters
    clusters(clusters==0) = []; % remove rejected cluster (indexed by zeros)
    nclust = numel(clusters);
    spike.label     = cell(1,nclust);
    spike.unit      = cell(1,nclust);
    spike.waveform  = cell(1,nclust);
    spike.timestamp = cell(1,nclust);
    ntt.TimeStamp = ntt.TimeStamp - ntt.hdr.FirstTimeStamp;
    for cl = 1:nclust
      unit_idx                  = ntt.CellNumber==clusters(cl);
      spike.label{cl}           = ['base-' num2str(cl)];
      spike.timestamp{cl}       = double(ntt.TimeStamp(unit_idx));
      spike.waveform{cl}        = ntt.dat(:,:,unit_idx);
      spike.unit{cl}            = clusters(cl).* ones(1,length(spike.timestamp{cl}));
    end
%     if isfield(ntt.hdr, 'NLX_Base_Class_Name')
%       spike.label   = {ntt.hdr.NLX_Base_Class_Name};
%     else
% %       spike.label   = {ntt.hdr.AcqEntName};
%         spike.label = {'base'}
%     end
%     spike.timestamp = {ntt.TimeStamp};
%     spike.waveform  = {ntt.dat};
%     spike.unit      = {ntt.CellNumber};
    spike.hdr       = ntt.hdr;
    fprintf('note ntt timestamps are in microseconds, first timestamp subtract so they start at 0.\n')

    
  case 'neuralynx_nse'
    % single channel file, read all records
    nse = read_neuralynx_nse(filename);
    if isfield(nse.hdr, 'NLX_Base_Class_Name')
      spike.label   = {nse.hdr.NLX_Base_Class_Name};
    else
      spike.label   = {nse.hdr.AcqEntName};
    end
    spike.timestamp = {nse.TimeStamp};
    spike.waveform  = {nse.dat};
    spike.unit      = {nse.CellNumber};
    spike.hdr       = nse.hdr;
    
  case 'neuralynx_nst'
    % single channel stereotrode file, read all records
    nst = read_neuralynx_nst(filename, 1, inf);
    if isfield(nst.hdr, 'NLX_Base_Class_Name')
      spike.label   = {nst.hdr.NLX_Base_Class_Name};
    else
      spike.label   = {nst.hdr.AcqEntName};
    end
    spike.timestamp = {nst.TimeStamp};
    spike.waveform  = {nst.dat};
    spike.unit      = {nst.CellNumber};
    spike.hdr       = nst.hdr;
    
  case 'neuralynx_nts'
    % single channel file, read all records
    nts = read_neuralynx_nts(filename);
    if isfield(nte.hdr, 'NLX_Base_Class_Name')
      spike.label   = {nts.hdr.NLX_Base_Class_Name};
    else
      spike.label   = {nts.hdr.AcqEntName};
    end
    spike.timestamp = {nts.TimeStamp(:)'};
    spike.waveform  = {zeros(0,length(nts.TimeStamp))};  % does not contain waveforms
    spike.unit      = {zeros(0,length(nts.TimeStamp))};  % does not contain units
    spike.hdr       = nts.hdr;
    
  case 'plexon_nex'
    % a single file can contain multiple channels of different types
    hdr  = read_plexon_nex(filename);
    typ  = [hdr.VarHeader.Type];
    chan = 0;
    
    spike.label     = {};
    spike.timestamp = {};
    spike.waveform  = {};
    spike.unit      = {};
    
    for i=1:length(typ)
      if typ(i)==0
        % neurons, only timestamps
        nex = read_plexon_nex(filename, 'channel', i);
        nspike = length(nex.ts);
        chan = chan + 1;
        spike.label{chan}     = deblank(hdr.VarHeader(i).Name);
        spike.waveform{chan}  = zeros(0, nspike);
        spike.unit{chan}      = nan(1,nspike);
        spike.timestamp{chan} = nex.ts;
      elseif typ(i)==3
        % neurons, timestamps and waveforms
        nex = read_plexon_nex(filename, 'channel', i);
        chan = chan + 1;
        nspike = length(nex.ts);
        spike.label{chan}     = deblank(hdr.VarHeader(i).Name);
        spike.waveform{chan}  = permute(nex.dat,[3 1 2]);
        spike.unit{chan}      = nan(1,nspike);
        spike.timestamp{chan} = nex.ts;
      end
    end
    spike.hdr = hdr;
    
  case 'plexon_nex5'
    % a single file can contain multiple channels of different types
    hdr  = read_nex5(filename);
    typ  = [hdr.VarHeader.Type];
    chan = 0;
    
    spike.label     = {};
    spike.timestamp = {};
    spike.waveform  = {};
    spike.unit      = {};
    
    for i=1:length(typ)
      if typ(i)==0
        % neurons, only timestamps
        nex = read_nex5(filename, 'channel', i);
        nspike = length(nex.ts);
        chan = chan + 1;
        spike.label{chan}     = deblank(hdr.VarHeader(i).Name);
        spike.waveform{chan}  = zeros(0, nspike);
        spike.unit{chan}      = nan(1,nspike);
        % spike.timestamp{chan} are the raw timestamps as recorded by the hardware system
        spike.timestamp{chan} = nex.ts;
      elseif typ(i)==3
        % waveform variables: timestamps and waveforms
        nex = read_nex5(filename, 'channel', i);
        chan = chan + 1;
        nspike = length(nex.ts);
        spike.label{chan}     = deblank(hdr.VarHeader(i).Name);
        spike.waveform{chan}  = permute(nex.dat,[3 1 2]);
        spike.unit{chan}      = nan(1,nspike);
        % spike.timestamp{chan} are the raw timestamps as recorded by the hardware system
        spike.timestamp{chan} = nex.ts;
      end
    end
    spike.hdr = hdr;
    
  case 'plexon_plx'
    % read the header information
    hdr   = read_plexon_plx(filename);
    nchan = length(hdr.ChannelHeader);
    typ   = [hdr.DataBlockHeader.Type];
    unit  = [hdr.DataBlockHeader.Unit];
    chan  = [hdr.DataBlockHeader.Channel];
    
    for i=1:nchan
      % select the data blocks that contain spike waveforms and that belong to this channel
      sel = (typ==1 & chan==hdr.ChannelHeader(i).Channel);
      
      if any(sel)
        % get the timestamps that correspond with this spike channel
        tsl = [hdr.DataBlockHeader(sel).TimeStamp];
        tsh = [hdr.DataBlockHeader(sel).UpperByteOf5ByteTimestamp];
        % convert the 16 bit high timestamp into a 32 bit integer
        ts = timestamp_plexon(tsl, tsh);
        spike.timestamp{i} = ts;
        spike.unit{i}      = unit(sel);
      else
        % this spike channel is empty
        spike.timestamp{i} = [];
        spike.unit{i}      = [];
      end
    end
    for i=1:nchan
      spike.label{i}    = deblank(hdr.ChannelHeader(i).Name);
      spike.waveform{i} = permute(read_plexon_plx(filename, 'ChannelIndex', i, 'header', hdr),[3 1 2]);
    end
    spike.hdr = hdr;
    
  case 'neuroshare' % NOTE: still under development
    % check that the required neuroshare toolbox is available
    ft_hastoolbox('neuroshare', 1);
    
    tmp = read_neuroshare(filename, 'readspike', 'yes');
    spike.label = {tmp.hdr.entityinfo(tmp.list.segment).EntityLabel};
    for i=1:length(spike.label)
      spike.waveform{i}  = tmp.spikew.data(:,:,i);
      spike.timestamp{i} = tmp.spikew.timestamp(:,i)';
      spike.unit{i}      = tmp.spikew.unitID(:,i)';
    end
    
  case 'neuroscope'
    % the information about the spikes is represented in:
    % x.clu.y or x.res.y (containing the timing +cluster info)
    % x.spk.y (containing the waveform info)
    % x.fet.y (containing features: do we need this?)
    
    if isfolder(filename)
      tmp = dir(filename);
      filenames = {tmp.name}';
    end
    
    % read the header
    filename_hdr = filenames{~cellfun('isempty',strfind(filenames,'.xml'))};
    hdr          = ft_read_header(fullfile(filename,filename_hdr), 'headerformat', 'neuroscope_xml');
    spikegroups  = hdr.orig.spikeGroups;
    fsample      = hdr.orig.rates.wideband;
    
    filename_clu = filenames(~cellfun('isempty',strfind(filenames,'.clu')));
    filename_spk = filenames(~cellfun('isempty',strfind(filenames,'.spk')));
    
    % FIXME should we do a sanity check on whether the clu and spk actually
    % belong together?
    
    c = cell(numel(filename_clu),1);
    w = cell(numel(filename_spk),1);
    for k = 1:numel(filename_clu)
      c{k} = LoadSpikeTimes(fullfile(filename,filename_clu{k}), fsample);
    end
    for k = 1:numel(filename_spk)
      w{k} = LoadSpikeWaveforms(fullfile(filename,filename_spk{k}),numel(spikegroups.groups{k}),spikegroups.nSamples(k));
    end
    
    spike = [];
    spike.label = cell(hdr.orig.spikeGroups.nGroups,1);
    spike.hdr   = hdr;
    spike.unit      = cell(1,numel(spike.label));
    spike.waveform  = cell(1,numel(spike.label));
    spike.timestamp = cell(1,numel(spike.label));
    
    for k = 1:numel(spike.label)
      sel = find(c{k}(:,3)>1); % values >1 corresponds to individual units, 0 = noise, 1 = MUA
      
      % the times are defined in s, convert to original time stamps
      timestamps = c{k}(sel,1) * hdr.orig.rates.wideband;
      if any(abs(timestamps-round(timestamps))>1e-5)
        ft_error('there seems to be a mismatch between the spike times and the expected integer-valued timestamps');
      end
      
      spike.timestamp{k} = round(timestamps(:))';
      spike.waveform{k}  = permute(w{k}(sel,:,:), [2 3 1]);
      spike.unit{k}      = c{k}(sel,3)';
      spike.label{k}     = sprintf('spikegroup%03d',k);
    end
    
  case 'nwb'
    ft_hastoolbox('MatNWB', 1);
    spike = read_nwb_spike(filename);
    
  case {'blackrock_nev'}
    % use the NPMK toolbox for the file reading
    ft_hastoolbox('NPMK', 1);
    
    % ensure that the filename contains a full path specification,
    % otherwise the low-level function fails
    [p, f, x] = fileparts(filename);
    if ~isempty(p)
      % this is OK
    elseif isempty(p)
      filename = which(filename);
    end
    
    % 'nosave' prevents the automatic conversion of the .nev file as a .mat file
    nev = openNEV(filename, 'nosave');
    
    nchan = length(nev.ElectrodesInfo);
    for i=1:nchan
      spike.label{i} = deblank(nev.ElectrodesInfo(i).ElectrodeLabel(:)');
      if isfield(nev.Data, 'Spikes')
        % select the spikes that were detected on this electrode
        sel = nev.Data.Spikes.Electrode;
        spike.timestamp{i} = nev.Data.Spikes.TimeStamp(sel);
        spike.unit{i}      = nev.Data.Spikes.Unit(sel);
        spike.waveform{i}  = nev.Data.Spikes.Waveform(sel);
      end
    end
    
  otherwise
    ft_error(['unsupported data format (' spikeformat ')']);
end

% add the waveform
if isfield(spike,'waveform')
  spike.dimord = '{chan}_lead_time_spike';
end
