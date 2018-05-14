classdef VariableMaskT4ResolveBuilder < mlfourdfp.CompositeT4ResolveBuilder
	%% VARIABLEMASKT4RESOLVEBUILDER  

	%  $Revision$
 	%  was created 13-May-2018 22:36:04 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
 		vmRBLog
 	end

	methods 
        
        %% OVERRIDES
		  
        function [ipr,this]   = revise(this, ipr)
            ipr = this.copySourceToDest(ipr); % crop/copy ipr.source to ipr.dest       
            this.imageRegLog = loggerFilename( ...
                ipr.dest{this.indexOfReference}, ...
                'func', 'CompositeT4ResolveBuilder_imageReg',...
                'path', ipr.logPath);
            this.resolveLog = loggerFilename( ...
                ipr.dest{this.indexOfReference}, ...
                'func', 'CompositeT4ResolveBuilder_t4ResolveAndPaste', ...
                'path', ipr.logPath);
            this.vmRBLog = loggerFilename( ...
                ipr.dest{this.indexOfReference}, ...
                'func', 'VariableMaskT4ResolveBuilder',...
                'path', ipr.logPath);
            
            this = this.imageReg(ipr);
            [ipr,~,this] = this.resolveAndPaste(ipr); 
            this.teardownRevision(ipr);
            this.rnumber = this.rnumber + 1;
        end  
        function this =         imageReg(this, ipr)
            stagedImgs  = this.lazyStageImages(ipr);    % contracted wrt this.indicesLogical
            blurredImgs = this.lazyBlurImages(ipr);     % "
            maskedImgs  = this.lazyMasksForImages(ipr); % "
            assertSizeEqual(stagedImgs, blurredImgs, maskedImgs);
            len = sum(this.indicesLogical);
            t4Failures = zeros(len, len);
            for m = 1:len
                for n = 1:len
                    if (m ~= n)
                        try
                            t4 = this.buildVisitor.filenameT4(stagedImgs{n}, stagedImgs{m});
                            if (~lexist(t4))
                                this.buildVisitor.align_crossModal( ...
                                    'dest',       blurredImgs{m}, ...
                                    'source',     blurredImgs{n}, ...
                                    'destMask',   'none', ...
                                    'sourceMask', 'none', ...
                                    't4',         t4, ...
                                    't4img_4dfp', false, ...
                                    'useMetricGradient', false, ...
                                    'log',        this.imageRegLog);
                                this.buildVisitor.align_crossModal( ...
                                    'dest',       blurredImgs{m}, ...
                                    'source',     blurredImgs{n}, ...
                                    'destMask',   maskedImgs{m}, ...
                                    'sourceMask', 'none', ...
                                    't4',         t4, ...
                                    't4img_4dfp', false, ...
                                    'useMetricGradient', true, ...
                                    'log',        this.imageRegLog);
                            end
                            % t4_resolve requires an idiomatic naming convention for t4 files,
                            % based on the names of image files
                            % e. g., image1_to_image2_t4    
                        catch ME
                            t4Failures(m,n) = t4Failures(m,n) + 1;
                            copyfile( ...
                                this.buildVisitor.transverse_t4, ...
                                this.buildVisitor.filenameT4(stagedImgs{n}, stagedImgs{m}), 'f');
                            dispwarning(ME);
                        end                        
                    end
                end
            end 
            this.assessT4Failures(t4Failures);
            this = this.assessT4ResolveErr(stagedImgs);
            this.deleteTrash;
        end
        
        %%
         
        function        assessT4Failures(this, t4fails)
            t4fails = sum(t4fails, 1);
            this.appendVMRBLog( ...
                sprintf('VariableMaskT4ResolveBuilder.assessT4Failures: t4Failures->%s\n', ...
                mat2str(t4fails)));
            this.indicesLogical(this.indicesLogical) = ...
                ensureRowVector(this.indicesLogical(this.indicesLogical)) & ...
                ensureRowVector(~logical(t4fails));   
            this.appendVMRBLog( ...
                sprintf('VariableMaskT4ResolveBuilder.assessT4Failures: this.indicesLogical->%s\n', ...
                mat2str(this.indicesLogical)));            
        end
        function this = assessT4ResolveErr(this, stagedImgs)
            len = sum(this.indicesLogical);
            this.t4_resolve_err = nan(len, len);
            for m = 1:length(stagedImgs)
                for n = 1:length(stagedImgs)
                    if (m ~= n)
                        try
                            this.t4_resolve_err(m,n) = ...
                                this.t4_resolve_errPairParser(mybasename(stagedImgs{m}), mybasename(stagedImgs{n}));
                        catch ME
                            dispwarning(ME);
                        end
                        
                    end
                end
            end 
            this.appendVMRBLog( ...
                sprintf('VariableMaskT4ResolveBuilder.assessT4ResolveErr: stagedImgs->%s\n', ...
                cell2str(stagedImgs, 'AsRow', true)));
            this.appendVMRBLog( ...
                sprintf('VariableMaskT4ResolveBuilder.assessT4ResolveErr: this.t4_resolve_err->%s\n', ...
                mat2str(this.t4_resolve_err)));            
        end
        function        appendVMRBLog(this, s)
            if (~lexist(this.vmRBLog))
                mlbash(sprintf('touch %s', this.vmRBLog));
            end
            mlbash(sprintf('echo \"%s\" >> %s', s, this.vmRBLog));
            fprintf(s);
        end        
        
 		function this = VariableMaskT4ResolveBuilder(varargin)
 			%% VARIABLEMASKT4RESOLVEBUILDER
 			%  @param .

 			this = this@mlfourdfp.CompositeT4ResolveBuilder(varargin{:});
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

