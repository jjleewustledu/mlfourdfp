classdef T4ResolveReporter 
	%% T4RESOLVEREPORTER 
    %  @param parameters needed by ctor.
    %  @return cell-array of instances of mlfourdfp.T4ResolveReport.  

	%  $Revision$
 	%  was created 04-May-2016 23:14:05
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.341360 (R2016a) for MACI64.
 	

	properties (Dependent)
 		reports
        sessionData
    end
    
    methods %% GET
        function r = get.reports(this)
            r = this.reports_;
        end
        function s = get.sessionData(this)
            s = this.sessionData_;
        end
    end

	methods 
		  
 		function this = T4ResolveReporter(varargin)
 			%% T4RESOLVEREPORTER
 			%  Usage:  this = T4ResolveReporter()

 			ip = inputParser;
            addParameter(ip, 'sessionData', this.sessionData_, @(x) isa(x, 'mlpipeline.ISessionData'));
            addParameter(ip, 'imagingPath', pwd,               @isdir);
            addParameter(ip, 'loggingPath', pwd,               @isdir);
            addParameter(ip, 'imagingFileprefix', '',          @ischar);
            addParameter(ip, 'loggingFileprefix', '',          @ischar);
            addParameter(ip, 'frameLength', 0,                 @isnumeric);
            parse(ip, varargin{:});
            
            import mlsystem.*;
            this.sessionData_ = ip.Results.sessionData;
            this.imagingDirTool_ = DirTool( ...
                fullfile(ip.Results.imagingPath, [ip.Results.imagingFileprefix '*.4dfp.img']));
            this.loggingDirTool_ = DirTool( ...
                fullfile(ip.Results.loggingPath, [ip.Results.loggingFileprefix '*.log']));
            
            import mlfourdfp.*;
            this.reports_ = cell(1, this.loggingDirTool_.length);
            for d = 1:this.loggingDirTool_.length
                this.reports_{d} = T4ResolveReport( ...
                    T4ResolveParser( ...
                        'sessionData',     this.sessionData, ...
                        'imagingFilename', '', ...
                        'loggingFilename', this.loggingDirTool_.fqfns{d}, ...
                        'frameLength', ip.Results.frameLength)); % this.strfindImagingFilename(this.loggingDirTool_.fqfns{d})
            end
 		end
    end 
    
    %% PROTECTED
    
    properties (Access = protected)
        sessionData_
        imagingDirTool_
        loggingDirTool_
        reports_
    end
    
    methods (Access = protected)
        function imgFn = strfindImagingFilename(this, logFn)
            for d = 1:this.imagingDirTool_.length
                [~,fp] = myfileparts(this.imagingDirTool_.fns{d});
                if (lstrfind(logFn, fp))
                    imgFn = this.imagingDirTool_.fqfns{d};
                    return
                end
            end
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

