classdef T4ResolveError < mlfourdfp.AbstractT4ResolveError
	%% T4RESOLVEERROR  

	%  $Revision$
 	%  was created 14-May-2018 15:13:51 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	    
    methods (Static)
        function em    = errorMat(varargin)
            this = mlfourdfp.T4ResolveError(varargin{:});
            [~,em] = this.estimateErr(this.theImages);
        end
    end
    
	methods
        
        %%
        
        function s     = summarizeFrames(this)
            for e = 1:this.sessionData.supEpoch
                this.sessionData.epoch = e;
                pwd0 = pushd(this.sessionData.tracerLocation);
                if (isempty(dir('*_t4'))); break; end
                [simgs,this] = this.stagedImgs(this.sessionData);
                [~,s{e}] = this.estimateErr(simgs); %#ok<AGROW>
                popd(pwd0);
            end
        end
        function [imgs,this] = stagedImgs(this, sd)
            %% STAGEDIMGS will not replace any existing indicesLogical
            
            ipr.dest = sd.tracerRevision('typ','fqfp');
            imgs = {};
            f1 = 1;
            this.indicesLogical = this.assessValidFrames(ipr);
            for f = 1:length(this.indicesLogical)
                if (this.indicesLogical(f))
                    ipr.currentIndex = f;
                    imgs{f1} = this.lazyExtractImage(ipr); %#ok<AGROW>
                    f1 = f1 + 1;
                end
            end
            assert(length(imgs) == sum(this.indicesLogical));
        end 
        function idx   = assessValidFrames(this, ipr)
            if (~lstrfind(ipr.dest, '.4dfp.ifh'))
                ipr.dest = [ipr.dest '.4dfp.ifh'];
            end
            d = mlfourd.ImagingContext(ipr.dest);
            if (isa(d, 'mlfourd.ImagingContext'))
                d = d.numericalNiftid;
            end
            if (isa(d, 'mlfourd.NIfTId'))
                d = mlfourd.NumericalNIfTId(d);
            end
            
            d   = d.volumeAveraged;
            idx = d.img > this.sessionData.fractionalImageFrameThresh * median(d.img) + this.noiseFloorOfCounts;
            idx = ensureRowVector(idx);
            this.logger.add('mlfourdfp.T4ResolveError.assessValidFrames.idx->%s\n', mat2str(idx));
        end
		  
 		function this = T4ResolveError(varargin)
 			%% T4RESOLVEERROR
 			%  @param .
            
 			this = this@mlfourdfp.AbstractT4ResolveError(varargin{:});
 		end
    end
    
    %% PROTECTED
    
    methods (Access = protected)
        function fp    = fileprefixIndexed(~, ipr)
            assert(isfield(ipr, 'dest'));
            assert(isfield(ipr, 'currentIndex'));
            fp = sprintf('%s_frame%i', ipr.dest, ipr.currentIndex);
        end
        function fqfp  = lazyExtractImage(this, ipr)
            %% LAZYEXTRACTIMAGE uses specifiers in ipr; will not replace any existing frame
            
            fqfp = this.fileprefixIndexed(ipr);
            if (~this.buildVisitor_.lexist_4dfp(fqfp))
                this.buildVisitor_.extract_frame_4dfp(ipr.dest, ipr.currentIndex, ['-o' fqfp]);
            end
        end
        function anImg = representativeImgs(this)
            if (iscell(this.theImages_))
                anImg = cell2str(this.theImages_, 'AsRow', true);
                if (length(anImg) > 79)
                    anImg = [anImg(1:79) '_'];
                end
                return
            end
            anImg = this.theImages_;
        end
        function this  = updateLogging(this)
            import mlpipeline.*;
            this = this.setLogPath(fullfile(pwd, 'Log', '')); % See also meanAndStd, summarize*
            ensuredir(this.getLogPath);
            this.logger_ = Logger( ...
                Logger.loggerFileprefix( ...
                    this.representativeImgs, ...
                    'func', 'T4ResolveErr',...
                    'path', this.getLogPath));
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

