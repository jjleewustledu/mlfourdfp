classdef CompositeT4ResolveError < mlfourdfp.AbstractT4ResolveError
	%% COMPOSITET4RESOLVEERROR  

	%  $Revision$
 	%  was created 29-May-2018 18:13:44 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
    methods (Static)
        function em    = errorMat(varargin)
            this = mlfourdfp.CompositeT4ResolveError(varargin{:});
            [~,em] = this.estimateErr(this.theImages);
        end
    end

	methods 
        
        %%
		  
        function s     = summarizeComposite(this)
            pwd0 = pushd(this.sessionData.tracerLocation);
            [simgs,this] = this.stagedImgs(this.sessionData);
            [~,s] = this.estimateErr(simgs);
            popd(pwd0);
        end
        function [imgs,this] = stagedImgs(this, varargin)
            ipr.dest = this.theImages;
            assert(iscell(ipr.dest));
            this = this.ensureSizeOfIndicesLogical(ipr.dest);
            imgs = {};
            d1 = 1;
            for d = 1:length(ipr.dest)
                if (this.indicesLogical(d))
                    imgs{d1} = ipr.dest{d}; %#ok<AGROW>
                    d1 = d1 + 1;
                end
            end
        end  
        
 		function this = CompositeT4ResolveError(varargin)
 			%% COMPOSITET4RESOLVEERROR
 			%  @param .

 			this = this@mlfourdfp.AbstractT4ResolveError(varargin{:});
 		end
    end 

    %% PROTECTED
    
    methods (Access = protected)
        function anImg = representativeImgs(this)
            assert(iscell(this.theImages_))
            anImg = '';
            for i = 1:length(this.theImages_)
                anImg = [anImg this.theImages_{i} '_']; %#ok<AGROW>
                if (i > 10); break; end
            end
            anImg = anImg(1:end-1);
        end
        function this  = updateLogging(this)
            import mlpipeline.*;
            this = this.setLogPath(fullfile(pwd, 'Log', '')); % See also meanAndStd, summarize*
            ensuredir(this.getLogPath);
            this.logger_ = Logger( ...
                fullfile(this.getLogPath, ...
                sprintf('%s_CompositeT4ResolveErr_%s', this.representativeImgs, mydatetimestr(now))));
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

