classdef T4ResolveError < mlfourdfp.AbstractT4ResolveError
	%% T4RESOLVEERROR  

	%  $Revision$
 	%  was created 14-May-2018 15:13:51 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	    
    methods (Static)
        function em = errorMat(varargin)
            this = mlfourdfp.T4ResolveError(varargin{:});
            [~,em] = this.estimateErr(this.theImages);
        end
        function emm = errorMatEpochs(varargin)
            this = mlfourdfp.T4ResolveError(varargin{:});
            sess = this.sessionData;
            Nt = length(sess.times);
            Ne = sess.maxLengthEpoch;
            emm = nan(Nt, Nt);
            pwd0 = pushd(sess.tracerLocation);
            for e = 1:sess.supEpoch-1
                pwde = pushd(fullfile(sprintf('E%i', e), 'Log', ''));
                emm((e-1)*Ne+1:e*Ne) = readBlockDiag(e);
                popd(pwde);
            end
            popd(pwd0);
            for c = 1:Nt
                emm(c,c) = 0;
            end
            csvwrite(emm, fullfile(sess.sessionPath, ...
                [sess.tracerRevision('typ','fp') '_T4ResolveErr_errorMatEpochs.csv']));
            
            function em = readBlockDiag(e_) 
                try      
                    dt = mlsystem.DirTool(sprintf('%se%ir*_T4ResolveErr_*.mat', lower(sess.tracer), e_));
                    [~,idx] = max(dt.itsListing.datenum);
                    load(dt.fns{idx}, 'em');              
                catch ME
                    dispexcept(ME, 'mlfourdfp:RuntimeError', ...
                        'T4ResolveError.errorMatEpochs.readBlockDiag had trouble examining %s', pwd);
                end
            end
        end
    end
    
	methods
        
        %%
        
        function idx         = assessValidFrames(this, ipr)
            if (~lstrfind(ipr.dest, '.4dfp.hdr'))
                ipr.dest = [ipr.dest '.4dfp.hdr'];
            end
            d   = mlfourd.ImagingContext2(ipr.dest);            
            d   = d.volumeSummed;
            n   = d.nifti;
            idx = n.img > this.sessionData.fractionalImageFrameThresh * median(n.img) + ...
                          mlnipet.NipetRegistry.instance().noiseFloorOfActivity;
            idx = ensureRowVector(idx);
            this.logger.add('mlfourdfp.T4ResolveError.assessValidFrames.idx->%s\n', mat2str(idx));
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
        function s           = summarizeFrames(this)
            for e = 1:this.sessionData.supEpoch
                this.sessionData.epoch = e;
                pwd0 = pushd(this.sessionData.tracerLocation);
                if (isempty(dir('*_t4'))); break; end
                [simgs,this] = this.stagedImgs(this.sessionData);
                [~,s{e}] = this.estimateErr(simgs); %#ok<AGROW>
                popd(pwd0);
            end
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
                anImg = strrep(cell2str(mybasename(this.theImages_), 'AsRow', true), ' ', '__');
                if (length(anImg) > 79)
                    anImg = [anImg(1:79) '_'];
                end
                return
            end
            anImg = this.theImages_;
        end
        function this  = updateLogging(this)
            import mlpipeline.*;
            [pth,fp] = myfileparts(this.representativeImgs);
            logpth = fullfile(pth, 'Log', '');
            ensuredir(logpth);            
            this.logger_ = Logger( ... ...
                fullfile(logpth, ...
                sprintf('%s_T4ResolveErr_%s', fp, mydatetimestr(now))));
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

