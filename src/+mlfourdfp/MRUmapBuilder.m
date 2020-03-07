classdef (Abstract) MRUmapBuilder < mlfourdfp.AbstractSessionBuilder & mlfourdfp.IUmapBuilder
	%% MRUMAPBUILDER  

	%  $Revision$
 	%  was created 22-Apr-2019 19:03:14 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.5.0.1067069 (R2018b) Update 4 for MACI64.  Copyright 2019 John Joowon Lee.
 	
	methods (Abstract)
 		s = mrSeriesLabel(this)
 	end

	methods 
        function     teardownLogs(this)
            ensuredir(this.getLogPath);
            try
                movefiles('*.log', this.getLogPath); 
                movefiles('*.txt', this.getLogPath);   
                movefiles('*.lst', this.getLogPath);    
                movefiles('*.mat0', this.getLogPath);   
                movefiles('*.sub', this.getLogPath); 
            catch ME
                dispwarning(ME, 'mlfourdfp:RuntimeWarning', ...
                    'PseudoCTBuilder.teardownLogs failed to move files into %s', this.getLogPath);
            end
        end
        function     teardownT4s(this)
            if (this.keepForensics); return; end
        end 
        function d = umapDicomPath(this)
            dt = mlsystem.DirTool( ...
                fullfile(this.sessionData.sessionPath, 'umaps', [this.mrSeriesLabel '_DT*.*'], ''));
            assert(dt.length > 0);
            d = dt.fqdns{1};
        end
        function d = umapPath(this)
            d = fileparts(this.umapDicomPath);
        end
        function umap = umapSynth(this, varargin)
            %% is this class' naming convention for mlpipeline.ISessionData.umapSynth.
            
            umap = this.sessionData.umapSynth('tracer', '', 'blurTag', '', 'typ', 'fqfp', varargin{:});
        end
		  
 		function this = MRUmapBuilder(varargin)
 			%% MRUMAPBUILDER
 			%  @param .

 			this = this@mlfourdfp.AbstractSessionBuilder(varargin{:});
            this.sessionData.attenuationCorrected = false;
            this.finished_ = mlpipeline.Finished(this, ...
                'path', this.getLogPath, 'tag', lower(this.sessionFolder));
 		end
    end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

