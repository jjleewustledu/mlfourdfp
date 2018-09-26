classdef DicomSorter < mlpipeline.DicomSorter
	%% DICOMSORTER  
    %  Top level:  sessionSort.

	%  $Revision$
 	%  was created 18-Sep-2016 15:18:38
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.341360 (R2016a) for MACI64.
 	
    methods (Static)
        function this  = Create(varargin)
            this = mlfourdfp.DicomSorter(varargin{:});
        end  
        function this  = CreateSorted(varargin)
            this = mlfourdfp.DicomSorter(varargin{:});
            this = this.sessionSort(varargin{:});
        end  
        function [s,r] = copyConverted(varargin)
            fv = mlfourdfp.FourdfpVisitor;
            [s,r] = fv.copy_4dfp(varargin{:});
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
        function tf    = lexistConverted(fqfp)
            fv = mlfourdfp.FourdfpVisitor;
            tf = fv.lexist_4dfp(fqfp);
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
                    this = this.sessionDcmConvert( ...
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
        
        %%
        
        function pth     = ctDestPath(this, str)
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
        function canonFp = dcm2imagingFormat(this, info, parentFqdn)
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
                this.studyData.seriesDicomAsterisk(parentFqdn), ...
                'base', canonFp, ...
                'options', opts);
            this.appendInfoFieldToIfh(info, 'KVP', canonFp);
            this.appendInfoFieldToIfh(info, 'RescaleIntercept', canonFp);
            this.appendInfoFieldToIfh(info, 'RescaleSlope', canonFp);
        end
        function pth     = destPath(this, str)
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
        function g       = getDcmConverter(~)
            g = 'dcm_to_4dfp';
        end
        
        function this = DicomSorter(varargin)
 			%% DICOMSORTER
            
            this = this@mlpipeline.DicomSorter(varargin{:});
            this.buildVisitor_ = mlfourdfp.FourdfpVisitor;
 		end
 	end 
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

