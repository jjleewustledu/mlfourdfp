classdef DicomSorter 
	%% DICOMSORTER  
    %  Top level:  sessionSort.

	%  $Revision$
 	%  was created 18-Sep-2016 15:18:38
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.341360 (R2016a) for MACI64.
 	
    
    properties
        cachedDcminfosFilename = 'DicomSorter_dcminfos_infos.mat';
        preferredInfoFields = {'SequenceName' 'SeriesDescription' 'ProtocolName' 'ImageType' 'SequenceName'}
        resetCached = false
        scansTags = {'scans' 'SCANS'}
    end
    
	properties (Dependent)
 		buildVisitor
        dicomExtension
        sessionData
        studyData
    end
    
    methods %% GET/SET
        function this = set.buildVisitor(this, v)
            assert(isa(v, 'mlfourdfp.FourdfpVisitor'));
            this.buildVisitor_ = v;
        end
        function v    = get.buildVisitor(this)
            v = this.buildVisitor_;
        end
        function e    = get.dicomExtension(this)
            e = this.studyData.dicomExtension;
        end
        function sd   = get.sessionData(this)
            sd = this.sessionData_;
        end
        function sd   = get.studyData(this)
            assert(~isempty(this.sessionData_));
            sd = this.sessionData_.studyData;
        end
    end

    methods (Static)
        function re    = ctFolderRegexp(str)
            %% CTFOLDERREGEXP
            %  @param string for regexp.
            %  @returns regexp result with name:  subjid.
            
            assert(ischar(str));
            re = regexp(str, '(?<subjid>^([A-Za-z]+\d+|[A-Za-z]+\d+_\d+|[A-Za-z]+\d+-\d+))(_|-)[A-Za-z_-]+\w*', 'names');
        end        
        function pth   = findRawdataSession(sessd)
            %% FINDRAWDATASESSION is a pseudo-inverse to DicomSorter.destPath
            %  @returns path to rawdata session.
            
            assert(isa(sessd, 'mlpipeline.SessionData'));
            
            import mlsystem.* mlfourdfp.*;
            dt    = DirTool(fullfile(sessd.rawdataDir, '*')); 
            res   = cellfun(@(x) DicomSorter.folderRegexp(x),          dt.dns, 'UniformOutput', false);
            keys  = cellfun(@(x) sprintf('%s_V%s', x.subjid, x.visit), res,    'UniformOutput', false);
            keys  = cellfun(@DicomSorter.fillEmpty, keys,                      'UniformOutput', false);
            m     = containers.Map(keys, dt.dns);            
            pth   = m(sprintf('%s_V%i', sessd.sessionLocation('typ', 'folder'), sessd.vnumber));
            pth   = fullfile(sessd.rawdataDir, pth, '');
        end
        function re    = folderRegexp(str)
            %% FOLDERREGEXP
            %  @param string for regexp.
            %  @returns regexp result with name:  subjid, visit.
            
            assert(ischar(str));
            re = regexp(str, '(?<subjid>^([A-Za-z]+\d+|[A-Za-z]+\d+_\d+|[A-Za-z]+\d+-\d+))(_|-)[A-Za-z_-]+(?<visit>\d)\w*', 'names');
        end
        function this  = sessions_to_4dfp(varargin)
            ip = inputParser;
            addParameter(ip, 'sessionFilter', @(x) iscell(x) || ischar(x));
            addParameter(ip, 'seriesFilter',  @(x) iscell(x) || ischar(x));
            addParameter(ip, 'sessionData', [], @(x) isa(x, 'mlpipeline.SessionData'));
            addParameter(ip, 'preferredName', 'unknown', @ischar);
            parse(ip, varargin{:});
            
            import mlfourdfp.* mlsystem.*;
            this = DicomSorter('sessionData', ip.Results.sessionData);
            dt = DirTools(ip.Results.sessionFilter);
            for idns = 1:length(dt.dns)
                try
                    cd(this.sessionData.rawdataDir);
                    srcPth = dt.dns{idns};
                    this = this.session_to_4dfp( ...
                        srcPth, ...
                        this.destPath(srcPth), ...
                        'sessionData', this.sessionData, ...
                        'seriesFilter', ip.Results.seriesFilter, ...
                        'preferredName', ip.Results.preferredName);
                catch ME
                    handwarning(ME);
                end
            end
        end        
        function this  = session_to_4dfp(varargin)
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'srcPath',            @isdir); % top-level folder for session raw data
            addOptional( ip, 'destPath', pwd,      @ischar);
            addParameter(ip, 'sessionData', [],      @(x) isa(x, 'mlpipeline.SessionData'));
            addParameter(ip, 'seriesFilter', {[]}, @(x) iscell(x) || ischar(x));
            addParameter(ip, 'preferredInfoFields', {'SeriesDescription'}, @iscell);
            addParameter(ip, 'preferredName', 'unknown', @ischar);
            parse(ip, varargin{:});            
            if (~isdir(ip.Results.destPath))
                mkdir(ip.Results.destPath);
            end
            
            import mlfourdfp.* mlsystem.* mlio.*;
            this = DicomSorter('sessionData', ip.Results.sessionData);
            this.preferredInfoFields = ip.Results.preferredInfoFields;
            pwd0 = pwd;
            cd(ip.Results.srcPath);
            [infos,fqdns] = this.findDcmInfos(ip.Results.seriesFilter);
            
            filteredNames = {};
            fv = mlfourdfp.FourdfpVisitor;
            for idns = 1:length(fqdns)
                try
                    canonFp = this.dcm_to_4dfp(infos{idns}, fqdns{idns});
                    if (~fv.lexist_4dfp(fullfile(ip.Results.destPath, canonFp)))
                        movefile([canonFp '.4dfp.*'], ip.Results.destPath);
                    end
                    canonFqfp = fullfile(ip.Results.destPath, canonFp);
                    filteredNames = [filteredNames canonFqfp]; %#ok<AGROW>
                catch ME
                    handwarning(ME);
                end
            end
            if (~isempty(filteredNames))
                this.linkPreferredSeries(filteredNames, ip.Results.preferredName);
            end
            
            cd(pwd0);
        end
        function this  = sessionSort(varargin)  
            %% SESSIONSORT sorts t1_mprage_sag, t2_spc_sag, TOF.
            
            ip = inputParser;
            addRequired( ip, 'srcPath',       @isdir); % top-level folder for session raw data
            addOptional( ip, 'destPath', pwd, @ischar);
            addParameter(ip, 'sessionData', [], @(x) isa(x, 'mlpipeline.SessionData'));
            parse(ip, varargin{:});

            import mlfourdfp.*;
            if (lstrfind(ip.Results.srcPath, 'CT')  || ...
                lstrfind(ip.Results.srcPath, 'CAL') || ...
                strcmpi( ip.Results.sessionData.modality, 'ct'))
                seriesList = {'AC_CT'};
                targetList = {'ct'};
            else
                seriesList = {'t1_mprage_sag' 't2_spc_sag' 'TOF'};
                targetList = {'mpr' 't2' 'tof'};
            end
            %seriesList = DicomSorter.SERIES_LIST;
            for s = 1:length(seriesList)
                try
                this = DicomSorter.session_to_4dfp( ...
                    ip.Results.srcPath, ip.Results.destPath, ...
                    'seriesFilter', seriesList{s}, ...
                    'sessionData', ip.Results.sessionData, ...
                    'preferredName', targetList{s});
                catch ME
                    dispwarning(ME);
                end
            end
        end
        function         linkPreferredSeries(filteredNames, preferredName)
            assert(iscell(filteredNames));
            assert(ischar(preferredName));
            
            import mlfourdfp.*;
            pwd0 = pwd;
            pth = myfileparts(filteredNames{1});
            if (~isempty(pth)); cd(pth); end
            filteredNames = DicomSorter.sortBySeriesNumber(filteredNames);
            if (~FourdfpVisitor.lexist_4dfp(preferredName))
                fv = FourdfpVisitor;
                fv.copy_4dfp(filteredNames{1}, preferredName);
            end
            cd(pwd0);
        end
        function names = sortBySeriesNumber(names)
            seriesNums = zeros(1,length(names));
            for n = 1:length(names)                
                re = regexp(names{n}, '\S+_series(?<series>\d+)\S*', 'names');
                seriesNums(n) = str2double(re.series);
            end
            tbl   = table(seriesNums', names', 'VariableNames', {'seriesNums' 'names'});
            tbl   = sortrows(tbl);
            names = tbl.names;
        end
    end
    
	methods
        function [s,r]         = appendInfoFieldToIfh(~, info, ifield, canonFp)
            if (~isfield(info, ifield))
                return
            end
            if (isnumeric(info.(ifield)))
                str = sprintf('%s := %g', ifield, info.(ifield));
            else
                str = sprintf('%s := %s', ifield, info.(ifield));
            end
            [s,r] = mlbash(sprintf('echo "%s" >> %s.4dfp.ifh', str, canonFp));
        end
        function n             = canonicalName(this, info)
            %% CANONICALNAME 
            %  @param dcminfo result, a struct.
            %  @returns <info.SeriesDescription>_series<info.SeriesNumber>.
            
            assert(isstruct(info) && ~isempty(info));
            n = sprintf('%s_series%g', info.SeriesDescription, info.SeriesNumber);
            n = this.scrubCanonicalName(n);
        end
        function pth           = ctDestPath(this, str)
            %% CTDESTPATH
            %  @param string for processing by regexp.
            %  @returns well-formed path for CT images as 4dfp.
            
            re  = this.ctFolderRegexp(str);
            if (~isempty(re))
                if (~isempty(re.subjid))
                    pth = fullfile(this.sessionData.subjectsDir, this.scrubSubjectID(re.subjid), '');
                    return
                end
            end
            pth = pwd;
        end
        function canonFp       = dcm_to_4dfp(this, info, parentFqdn)
            %% DCM_TO_4DFP
            %  @param info is a struct produced by dcminfo.
            %  @param parentFqdn is the parent directory of a DICOM directory.
            %  @returns canonFp, a fileprefix determined by this.canonicalName.
            
            assert(isstruct(info) && ~isempty(info))
            assert(isdir(parentFqdn));

            canonFp = this.canonicalName(info);
            opts = '-g';
            if (lstrfind(canonFp, 'TOF'))
                opts = [opts ' -t T'];
            end
            this.buildVisitor_.dcm_to_4dfp( ...
                this.sessionData.seriesDicomAsterisk(parentFqdn), ...
                'base', canonFp, ...
                'options', opts);
            this.appendInfoFieldToIfh(info, 'KVP', canonFp);
            this.appendInfoFieldToIfh(info, 'RescaleIntercept', canonFp);
            this.appendInfoFieldToIfh(info, 'RescaleSlope', canonFp);
        end
        function [infos,fqdns] = dcmInfos(this, varargin)
            %% DCMINFOS
            %  @param srcPath points to rawdata/SessionName has the pwd as default.
            %  @returns infos, a cell array containing struct results of dcminfo acting in srcPath.
            %  @returns fqdns, a cell array containing /path/to/DICOM.
            
            ip = inputParser;
            addOptional(ip, 'srcPath', pwd, @isdir); % top-level folder, not 'SCANS'
            parse(ip, varargin{:});
            
            import mlsystem.* mlio.*;
            cd(ip.Results.srcPath);
            scans = DirTool('*');
            for iscan = 1:length(scans.fqdns)
                if (lstrfind(scans.fqdns{iscan}, this.scansTags))
                    cd(scans.fqdns{iscan});  
                    series = DirTool('*');
                    fqdns  = series.fqdns;
                    if (lexist(this.cachedDcminfosFilename, 'file') && ~this.resetCached)
                        load(this.cachedDcminfosFilename);
                        return
                    end                    
                    infos  = cell(1, length(fqdns));
                    for iseries = 1:length(fqdns)
                        try
                            dcms = DirTool(fullfile(fqdns{iseries}, 'DICOM', ['*.' this.dicomExtension]));
                            if (~isempty(dcms.fqfns))
                                infos{iseries} = dicominfo(dcms.fqfns{1});
                            end
                            save(this.cachedDcminfosFilename, 'infos');
                        catch ME
                            handwarning(ME);
                        end
                    end 
                end                
            end
        end
        function pth           = destPath(this, str)
            %% DESTPATH 
            %  @param str, a string for rawdata source folder.
            %  @returns pth, the path to a canonical session path determined by internally stored mlpipeline.SessionData
            %  object.
            %  Uses:  this.ctDestPath, this.folderRegexp, this.scrubSubjectID.
            
            assert(ischar(str));
            if (lstrfind(lower(str), 'ct'))
                pth = this.ctDestPath(str);
                return
            end
            
            re  = this.folderRegexp(str);
            if (~isempty(re))
                if (~isempty(re.subjid))
                    if (~isempty(re.visit))
                        pth = fullfile(this.sessionData.subjectsDir, this.scrubSubjectID(re.subjid), sprintf('V%s', re.visit), '');
                        return
                    else
                        pth = fullfile(this.sessionData.subjectsDir, this.scrubSubjectID(re.subjid), '');
                        return
                    end
                end
            end
            pth = pwd;
        end
        function tf            = filterMatch(this, varargin)
            %% FILTERMATCH
            %  @param info is a struct produced by dcminfo.
            %  @param tomatch is a filtering string.
            %  @param named 'fields' is a cell array of field names to inquire of info; default is
            %  this.preferredInfoFields.
            %  @returns tf, which is logical.
            
            ip = inputParser;
            addRequired( ip, 'info', @isstruct);
            addRequired( ip, 'tomatch', @(x) ~isempty(x) && ischar(x));
            addParameter(ip, 'fields', this.preferredInfoFields, @iscell);
            parse(ip, varargin{:});
            
            values = {};            
            for f = 1:length(ip.Results.fields)
                if (isfield(ip.Results.info, ip.Results.fields{f}))
                    value  = ip.Results.info.(ip.Results.fields{f});
%                     if (isstruct(value))
%                         value = struct2str(value);
%                     end
%                     if (isnumeric(value))
%                         value = num2str(value);
%                     end
                    if (ischar(value))
                        values = [values this.scrubCanonicalName(value)]; %#ok<AGROW>
                    end
                end
            end
            tf = lstrfind(values, ip.Results.tomatch);
        end
        function [info, fqdn ] = findDcmInfo(this, ordinal, tofind, varargin)
            %% FINDDCMINFO
            %  @param ordinal is the index of desireable results from findDcmInfos:  1, 2, 3, ....
            %  @param tofind is a string or cell-array of strings for filtering information from dcminfo operating in srcPath.
            %  @param srcPath points to rawdata/SessionName but has the pwd as default.
            %  @returns info, the struct obtained from dcminfo.
            %  @returns fqdn, the /path/to/DICOM.
            
            assert(isnumeric(ordinal));
            [infos,fqdns] = this.findDcmInfos(tofind, varargin{:});
            assert(~isempty(infos));
            assert(~isempty(fqdns));
            info = infos{ordinal};
            fqdn = fqdns{ordinal};
        end
        function [infos,fqdns] = findDcmInfos(this, tofind, varargin)
            %% FINDDCMINFOS            
            %  @param tofind is a string or cell-array of strings for filtering information from dcminfo operating in srcPath.
            %  @param srcPath points to rawdata/SessionName but has the pwd as default.
            %  @returns infos, a cell array of info structs obtained from dcminfo.
            %  @returns fqdns, a cell array containing /path/to/DICOM.
            
            if (ischar(tofind))
                tofind = {tofind};
            end
            assert(iscell(tofind));
            
            [infos0,fqdns0] = this.dcmInfos(varargin{:});
            infos = {};
            fqdns = {};
            for idx = 1:length(infos0)
                for fdx = 1:length(tofind)
                    if (~isempty(infos0{idx}))
                        if (this.filterMatch(infos0{idx}, tofind{fdx})) %%%, 'fields', fields(infos0{idx})))
                            infos = [infos infos0{idx}]; %#ok<AGROW>
                            fqdns = [fqdns fqdns0{idx}]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
        function n             = scrubCanonicalName(~, n)
            %% SCRUBCANONICALNAME removes '*', '[', ']', ' '; and replaces '<', '>', '(', ')', ':' with '_'.
            
            n = strrep(n, '*', '');
            n = strrep(n, '<', '_');
            n = strrep(n, '>', '_');
            n = strrep(n, '[', '');
            n = strrep(n, ']', '');
            n = strrep(n, '(', '_');
            n = strrep(n, ')', '_');
            n = strrep(n, ' ', '');
            n = strrep(n, ':', '_');
        end
        function n             = scrubSubjectID(~, n)
            %% SCRUBSUBJECTID replaces '-' with '_'.
            
            n = strrep(n, '-', '_');
        end
        
        function this = DicomSorter(varargin)
 			%% DICOMSORTER
 			%  Usage:  this = DicomSorter(param_name, param_value)
            %  @param 'studyData' is an mlpipeline.StudyDataHandle.
            %  left unset, DicomSorter attempts to use only sessionData.
            %  @param 'sessionData' is an mlpipeline.SessionData; 
            
            ip = inputParser;
            addParameter(ip, 'sessionData', [], @(x) isa(x, 'mlpipeline.SessionData')     || isempty(x));
            addParameter(ip, 'studyData',   [], @(x) isa(x, 'mlpipeline.StudyDataHandle') || isempty(x));
            parse(ip, varargin{:});
            
            this.sessionData_  = ip.Results.sessionData;
            if (~isempty(ip.Results.studyData))                
                this.sessionData_.studyData = ip.Results.studyData;
            end
            this.buildVisitor_ = mlfourdfp.FourdfpVisitor;
 		end
 	end 
    
    %% PRIVATE
    
    properties (Access = private)
        buildVisitor_
        sessionData_
    end
    
    methods (Static, Access = private)        
        function c = fillEmpty(c)
            if (isempty(c))
                c = num2str(rand);
            end
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

