classdef RawDataSorter 
	%% RAWDATASORTER:  start with dcm_sort_PPG; then use move* or copy*.

	%  $Revision$
 	%  was created 06-Sep-2016 15:36:57
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.341360 (R2016a) for MACI64.
 	
    
    properties
        dicomIndexUTE
    end
    
    properties (Dependent)
        studyData
        sessionData
        sessionPath
    end
    
    methods %% GET/SET
        function g    = get.studyData(this)
            assert(~isempty(this.sessionData_));
            g = this.sessionData.studyData;
        end
        function this = set.studyData(this, s)
            assert(isa(s, 'mlpipeline.StudyDataHandle'));
            assert(~isempty(this.sessionData_));
            this.sessionData_.studyData = s;
        end
        function g    = get.sessionData(this)
            g = this.sessionData_;
        end
        function this = set.sessionData(this, s)
            assert(isa(s, 'mlpipeline.SessionData'));
            this.sessionData_ = s;
        end
        function g    = get.sessionPath(this)
            g = this.sessionData_.sessionPath;
        end
        function this = set.sessionPath(this, s)
            assert(isdir(s));            
            this.sessionData_.sessionPath = s;
        end
    end
    
    methods (Static)
        function [s,r] = bash_copyfile(loc0, loc)
            if (isdir(loc))
                error('mlfourdfp:potentialDataCorruption', ...
                      'destination already exists:  RawDataSorter.bash_copyfile.loc->%s', loc);
            end
            if (~isdir(fileparts(loc)))
                mkdir(fileparts(loc));
            end
            try
                [s,r] = mlbash(sprintf('cp -r %s %s', loc0, fileparts(loc)));
                [~,fold0] = fileparts(loc0);
                loc1  = fullfile(fileparts(loc), fold0);
                if (~strcmp(loc1, loc))
                    [s,r] = mlbash(sprintf('mv %s %s', loc1, loc));
                end
            catch ME
                handexcept(ME);
            end
        end
        function [s,r] = bash_movefile(loc0, loc)
            if (isdir(loc))
                error('mlfourdfp:potentialDataCorruption', ...
                      'destination already exists:  RawDataSorter.bash_movefile.loc->%s', loc);
            end            
            if (~isdir(fileparts(loc)))
                mkdir(fileparts(loc));
            end
            try
                [s,r] = mlbash(sprintf('mv %s %s', loc0, fileparts(loc)));
                [~,fold0] = fileparts(loc0);
                loc1  = fullfile(fileparts(loc), fold0);
                if (~strcmp(loc1, loc))
                    [s,r] = mlbash(sprintf('mv %s %s', loc1, loc));
                end
            catch ME
                handexcept(ME);
            end
        end
        function this = copy(varargin)
            %% COPY is a convenience static function for copyRawData, copyUTE
            %  @params rawdataFolder from /path/to/rawdata.
            
            import mlfourdfp.* mlraichle.*;            
            ip = inputParser;
            addRequired( ip, 'rawdataFolder', @(x) isdir(RawDataSorter.sourceRawDataPath(x)) && ...
                                                   isdir(RawDataSorter.sourceScansPath(x)));
            parse(ip, varargin{:});
                       
            srcRawData  = RawDataSorter.sourceRawDataPath(ip.Results.rawdataFolder);
            iSrcFold    = strfind(srcRawData, 'rawdata/');
            iSrcFold    = iSrcFold + 8;
            re          = DicomSorter.folderRegexp(srcRawData(iSrcFold:end));
            studyd      = StudyData;
            sessd       = mlraichle.SessionData('studyData', studyd, 'sessionPath', RawDataSorter.destSessionPath(re.subjid));
            sessd.vnumber = str2double(re.visit);
            this        = RawDataSorter('sessionData', sessd);
            this        = this.dcm_sort_PPG(srcRawData);
            this        = this.copyRawData(srcRawData);
            this        = this.copyUTE( ...
                          RawDataSorter.sourceScansPath(ip.Results.rawdataFolder));
        end
        function this = copyConverted(varargin)
            
            import mlfourdfp.* mlraichle.*;
            
            ip = inputParser;
            addParameter(ip, 'session', '', @ischar); % e.g., $PPG/jjlee/HYGLY15, HYGLY15
            parse(ip, varargin{:});
            
            studyData = StudyData;
            sessPath  = RawDataSorter.destSessionPath(ip.Results.session);
            sessData  = SessionData('studyData', studyData, 'sessionPath', sessPath);
 			this      = RawDataSorter('sessionData', sessData);
            this      = this.copyRawDataConverted('sessionData', sessData, 'sessionPath', sessPath);
        end
        function this = move(varargin)
            %% MOVE is a convenience static function for moveRawData, copyUTE
            %  @params rawdataFolder from /path/to/rawdata.
            
            import mlfourdfp.* mlraichle.*;            
            ip = inputParser;
            addRequired( ip, 'oper',          @(x) isa(x, 'function_handle'));
            addRequired( ip, 'rawdataFolder', @(x) isdir(RawDataSorter.sourceRawDataPath(x)) && ...
                                                   isdir(RawDataSorter.sourceScansPath(x)));
            parse(ip, varargin{:});
                       
            srcRawData  = RawDataSorter.sourceRawDataPath(ip.Results.rawdataFolder);
            [~,srcFold] = fileparts(srcRawData);            
            re          = DicomSorter.folderRegexp(srcFold);         
            sessd       = SessionData('sessionPath', ...
                          RawDataSorter.destSessionPath(re.subjid));
            this        = RawDataSorter('sessionData', sessd);
            this        = this.moveRawData(srcRawData);
            this        = this.copyUTE( ...
                          RawDataSorter.sourceRawDataPath(ip.Results.rawdataFolder));
        end
        function pth  = sourceRawDataPath(str)
            assert(ischar(str));
            pth = fullfile(getenv('PPG'), 'rawdata', str, 'RESOURCES', 'RawData', '');            
            if(~isdir(pth))
                pth = ''; 
            end
        end
        function pth  = sourceScansPath(str)
            assert(ischar(str));
            if (~strcmp(str(1), '/'))
                pth = fullfile(getenv('PPG'), 'rawdata', str, 'SCANS', '');
                if (~isdir(pth))
                    pth = fullfile(getenv('PPG'), 'rawdata', str, 'scans', '');
                end
            else
                pth = fullfile(str, '..', '..', 'SCANS', '');
                if (~isdir(pth))
                    pth = fullfile(str, '..', '..', 'scans', '');
                end
            end
            if(~isdir(pth))
                pth = ''; 
            end
        end
        function pth  = destSessionPath(str)
            assert(ischar(str));
            if (~strcmp(str(1), '/'))
                pth = fullfile(mlraichle.RaichleRegistry.subjectsDir, str, '');
                assert(isdir(pth));
            else
                pth = str;
            end
            if(~isdir(pth))
                pth = ''; 
            end
        end
        function v    = vfind(str)
            v = 1;
            str = lower(str);
            strings = {'_v_' '_v' '_visit_' '_visit' '-v-'};
            for sidx = 1:length(strings)
                if (lstrfind(str, strings{sidx}))
                    idx = strfind(str, strings{sidx});
                    v   = str2double(str(idx + length(strings{sidx})));
                    if (~isnan(v) && v > 0)
                        break
                    end
                end
            end
            assert(~isempty(v));
            assert( isnumeric(v));
        end
    end
    
	methods 
        function this = copyUTE(this, varargin)
            %% COPYUTE
            %  @param srcLoc is the source directory, e.g., '/rawdata/HYGLY25_VISIT_1/SCANS'.
            %  @returns this.
            
            this = this.operationOnData(@mlfourdfp.RawDataSorter.bash_copyfile, @this.UTEMatch, varargin{:});
        end
        function this = copyRawData(this, varargin)
            %% COPYRAWDATA
            %  @param srcLoc is the source directory, e.g., '/rawdata/HYGLY25_VISIT_1/RESOURCES/RawData'.
            %  @returns this.
            
            this = this.operationOnData(@mlfourdfp.RawDataSorter.bash_copyfile, @this.rawDataMatch, varargin{:});
        end
        function this = copyRawDataConverted(this, varargin)
            %% COPYRAWDATACONVERTED copies legacy converted data from Lars Couture
            
            ip = inputParser;
            addParameter(ip, 'sessionData', this.sessionData, @(x) isa(x, 'mlpipeline.SessionData'));
            parse(ip, varargin{:});
                     
            this.sessionData_ = ip.Results.sessionData;
            [~,sfold] = fileparts(this.sessionData.sessionPath);

            for iv = 1:2
                vfold = sprintf('V%i', iv);
                vpath = fullfile(this.sessionData.sessionPath, vfold, '');
                cpath = fullfile(getenv('PPG'), 'converted', sfold, vfold, '');
                tracers = {'FDG_V' 'HO1_V' 'HO2_V' 'OC1_V' 'OC2_V' 'OO1_V' 'OO2_V'};
                for it = 1:length(tracers)
                    tfold = sprintf('%s%i', tracers{it}, iv);
                    oldName = fullfile(cpath, tfold, '');
                    newName = fullfile(vpath, tfold, '');
                    if (~isempty(newName))
                        this.binaryFilesystemOperation0(oldName, newName, @copyfile);
                    end
                end
            end
        end
        function this = dcm_sort(this, rawDataDir)
            assert(isdir(rawDataDir));
            dt = mlsystem.DirTool(rawDataDir);
            if (~isempty(dt.dns))
                return
            end
            vtor = mlfourdfp.FourdfpVisitor;
            vtor.dcm_sort(rawDataDir);
        end
        function this = dcm_sort_PPG(this, rawDataDir)
            assert(isdir(rawDataDir));
            dt = mlsystem.DirTool(rawDataDir);
            if (~isempty(dt.dns))
                return
            end
            vtor = mlfourdfp.FourdfpVisitor; % escalate to ad hoc shell script
            vtor.dcm_sort_PPG(rawDataDir);
        end
        function this = moveRawData(this, varargin)
            this = this.operationOnData(@mlfourdfp.RawDataSorter.bash_movefile, @this.rawDataMatch, varargin{:});
        end
        function [destLocs,srcLocs] = rawDataMatch(~, loc0, loc1)
            %% RAWDATAMATCH
            %  @param loc0 is the path created by dcm_sort for PET raw data, e.g., 
            %  '/rawdata/HYGLY25_VISIT_1/RESOURCES/RawData'
            %  @param loc1 is the path ending in the session visit folder, e.g., '/path/to/V1', '/path/to/V2', ....
            %  @returns destLocs, cell array of destination paths, e.g., {'/path/to/V1/FDG_V1' ...}.
            %  @returns loc0.
            
            assert(isdir(loc0));          
            [~,vfold] = fileparts(loc1);
            assert(lstrfind(vfold, 'V') && 2 == length(vfold));
            
            if (~isdir(loc1))
                [s,m,mid] = mkdir(loc1);
                if (~s)
                    error('mlfourdfp:mkdirFailed', 'RawDataSorter.rawDataMatch:  %s; %s', m, mid);
                end
            end
            dt = mlsystem.DirTool(loc0);
            srcLocs = dt.fqdns;
            destLocs = cell(size(srcLocs));
            for id = 1:length(srcLocs)
                fold = srcLocs{id};
                if (lstrfind(fold, 'Head_MRAC_PET_60min'))
                    destLocs{id} = fullfile(loc1, ['FDG_' vfold]);
                end
                if (lstrfind(fold, 'Head_MRAC_PET_5min'))
                    destLocs{id} = fullfile(loc1, ['Twilite_' vfold]);
                end
                if (lstrfind(fold, 'Head_HO'))
                    idxS = strfind(fold, 'HO');
                    destLocs{id} = fullfile(loc1, ['HO' fold(idxS+2) '_' vfold]);
                end
                if (lstrfind(fold, 'Head_OC'))
                    idxS = strfind(fold, 'OC');
                    destLocs{id} = fullfile(loc1, ['OC' fold(idxS+2) '_' vfold]);
                end
                if (lstrfind(fold, 'Head_OO'))
                    idxS = strfind(fold, 'OO');
                    destLocs{id} = fullfile(loc1, ['OO' fold(idxS+2) '_' vfold]);
                end
            end
        end
        function [destLocs,srcLocs] = UTEMatch(this, loc0, loc1)
            %% UTEMATCH
            %  @param loc0 is the path containing rawdata dicoms, e.g., '/rawdata/HYGLY25_VISIT_1/SCANS'.
            %  @param loc1 is the path ending in the session visit folder, e.g., '/path/to/V1', '/path/to/V2', ...
            %  @returns locs, a cell-array of locations, e.g., '/path/to/V1/FDG_V1/umap'.
            %  @returns fqdn, the path to source dicoms, e.g., '/rawdata/HYGLY25_VISIT_1/SCANS/45/DICOM'.

            tags = {'FDG_' 'HO1_' 'HO2_' 'OC1_' 'OC2_' 'OO1_' 'OO2_'};
            [destLocs0,srcLocs0] = umapMatch(this, loc0, loc1, 'MRAC_PET_5min_in_UMAP',  {'Twilite'});
            [destLocs1,srcLocs1] = umapMatch(this, loc0, loc1, 'UTE_AC_only_UMAP',       tags);
            [destLocs2,srcLocs2] = umapMatch(this, loc0, loc1, 'UTE_MRAC_UMAP',          tags);
            [destLocs3,srcLocs3] = umapMatch(this, loc0, loc1, 'MRAC_PET_60min_in_UMAP', tags);            
           
            destLocs = [destLocs0 destLocs1 destLocs2 destLocs3];
            srcLocs  = [srcLocs0  srcLocs1  srcLocs2  srcLocs3];
            
        end
        function [destLocs,srcLocs] = umapMatch(this, loc0, loc1, dcmInfoTag, tags)
            
            assert(isdir(loc0));          
            [~,vfold] = fileparts(loc1);
            assert(lstrfind(vfold, 'V') && 2 == length(vfold));            
            if (~isdir(loc1))
                [s,m,mid] = mkdir(loc1);
                if (~s)
                    error('mlfourdfp:mkdirFailed', 'RawDataSorter.dcminfoMatch:  %s; %s', m, mid);
                end
            end            
            ds = mlfourdfp.DicomSorter('sessionData', this.sessionData);
            try
                [info,fqdn] = ds.findDcmInfo(1, dcmInfoTag, fileparts(loc0));
            catch ME %#ok<NASGU>
                try                    
                    [info,fqdn] = ds.findDcmInfo(1, dcmInfoTag, fullfile(this.sessionData.rawdataDir, fileparts(loc0), ''));
                catch ME1 %#ok<NASGU>
                    fprintf('RawDataSorter.umapMatch:  no match for %s\n', dcmInfoTag);
                    info = [];
                end
            end
            srcLocs = {}; destLocs = {};
            if (~isempty(info))
                srcLocs  = cellfun(@(x) fullfile(fqdn, 'DICOM', ''),           tags, 'UniformOutput', false);
                destLocs = cellfun(@(x) fullfile(loc1, [x vfold], 'umap', ''), tags, 'UniformOutput', false);
            end
        end
		  
 		function this = RawDataSorter(varargin)
 			%% RAWDATASORTER
 			%  Usage:  this = RawDataSorter(param_name, param_value)
            %  @param 'studyData' is an mlpipeline.StudyDataHandle.
            %  @param 'sessionData' is an mlpipeline.SessionData.
            
            ip = inputParser;
            addParameter(ip, 'sessionData', [], @(x) isa(x, 'mlpipeline.SessionData'));
            addParameter(ip, 'studyData',   [], @(x) isa(x, 'mlpipeline.StudyDataHandle') || isempty(x));
            parse(ip, varargin{:});
            
            this.sessionData_ = ip.Results.sessionData;
            if (~isempty(ip.Results.studyData))
                this.sessionData_.studyData = ip.Results.studyData;
            end
 		end
    end 
    
    %% PRIVATE
    
    properties (Access = private)
        sessionData_
    end

    methods (Access = private)        
        function [s,m,mid] = binaryFilesystemOperation0(~, loc0, loc, funch)
            try
                assert(isdir(loc0));
                if (~lexist(fileparts(loc)))
                    [s,m,mid] = mkdir(fileparts(loc));
                    if (~s)
                        error('mlfourdfp:mkdirFailed', 'RawDataSorter.binaryFilesystemOperation:  %s; %s', m, mid);
                    end
                end
                if (~isdir(loc))
                    [s,m,mid] = funch(loc0, loc);
                    if (s)
                        fprintf('RawDataSorter.binaryFilesystemOperation:  transformed %s to %s\n', loc0, loc);
                    else
                        error('mlfourdfp:movefileFailed', 'RawDataSorter.binaryFilesystemOperation:  %s; %s', m, mid);
                    end
                end
            catch ME
                handwarning(ME);
            end
        end
        function [s,r] = binaryFilesystemOperation(~, loc0, loc, funch)
            try
                assert(isdir(loc0));
                [s,r] = funch(loc0, loc);
                fprintf('RawDataSorter.binaryFilesystemOperation:  transformed %s to %s\n', loc0, loc);
            catch ME
                handwarning(ME, 'RawDataSorter.binaryFilesystemOperation:  %s\n', r);
            end
        end
        function this = operationOnData(this, varargin)
            %% OPERATIONONDATA
            %  @param oper is a function_handle to the filesystem operation, e.g., copyfile, movefile.
            %  @param match is a function_handle to the matching function.
            %  @param srcLoc is the source directory interpreted by the dereference of match.
            %  @returns this after applying oper to results of [destLocs,srcLoc] = match.
            
            ip = inputParser;
            addRequired( ip, 'oper',  @(x) isa(x, 'function_handle'))
            addRequired( ip, 'match', @(x) isa(x, 'function_handle'))
            addRequired( ip, 'srcLoc', @isdir);
            parse(ip, varargin{:});

            assert(isa(this.sessionData, 'mlpipeline.SessionData'));
            [destLocs,srcLocs] = ip.Results.match( ...
                ip.Results.srcLoc, ...
                fullfile(this.sessionData.sessionPath, sprintf('V%i', this.sessionData.vnumber), ''));
            if (~iscell(srcLocs))
                srcLocs = {srcLocs}; 
            end
            if (~iscell(destLocs))
                destLocs = {destLocs}; 
            end
            assert(length(srcLocs) == length(destLocs));
            for d = 1:length(destLocs)
                try
                    this.binaryFilesystemOperation(srcLocs{d}, destLocs{d}, ip.Results.oper);
                catch ME
                    handwarning(ME);
                end
            end
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

