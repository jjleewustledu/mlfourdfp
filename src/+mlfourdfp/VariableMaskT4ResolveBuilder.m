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
            ipr.logPath = fullfile(this.sessionData.vallLocation, 'Log', '');
            ensuredir(ipr.logPath);
            ipr = this.copySourceToDest(ipr); % crop/copy ipr.source to ipr.dest       
            this.imageRegLog = this.loggerFilename( ...
                ipr.dest{this.indexOfReference}, ...
                'func', 'CompositeT4ResolveBuilder_imageReg',...
                'path', ipr.logPath);
            this.resolveLog = this.loggerFilename( ...
                ipr.dest{this.indexOfReference}, ...
                'func', 'CompositeT4ResolveBuilder_t4ResolveAndPaste', ...
                'path', ipr.logPath);
            this.vmRBLog = this.loggerFilename( ...
                ipr.dest{this.indexOfReference}, ...
                'func', 'VariableMaskT4ResolveBuilder',...
                'path', ipr.logPath);
            
            this.purgeT4s(ipr);
            
            ipr.cacheLoc = ensuredir( ...
                fullfile(this.sessionData.vallLocation, ['VarMaskT4RBCache_' mydatetimestr(now)]));
            ipr.maskForImages = 'none';
            ipr.useMetricGradient = true;
            this = this.imageReg(ipr);  
            [ipr,~,this] = this.resolveAndPaste(ipr); 
            
            ipr.cacheLoc = ensuredir( ...
                fullfile(this.sessionData.vallLocation, ['VarMaskT4RBCache_' mydatetimestr(now)]));
            ipr.maskForImages = 'Msktgen';
            ipr.useMetricGradient = true;
            this = this.imageReg(ipr);  
            [ipr,~,this] = this.resolveAndPaste(ipr); 
            
            ipr.cacheLoc = ensuredir( ...
                fullfile(this.sessionData.vallLocation, ['VarMaskT4RBCache_' mydatetimestr(now)]));
            ipr.maskForImages = 'Msktgen';
            ipr.useMetricGradient = true;
            this = this.imageReg(ipr);            
            [ipr,~,this] = this.resolveAndPaste(ipr); 
            
            this.teardownRevision(ipr);
            this.rnumber = this.rnumber + 1;
        end  
        function this =         imageReg(this, ipr)
            stagedImgs  = this.lazyStageImages(ipr);    % contracted wrt this.indicesLogical
            blurredImgs = this.lazyBlurImages(ipr);     % "
            maskedImgs  = this.lazyMasksForImages(ipr, stagedImgs); % "
            assertSizeEqual(stagedImgs, blurredImgs, maskedImgs);
            len = sum(this.indicesLogical);
            t4Failures = zeros(len, len);
            for m = 1:len
                for n = 1:len
                    if (m ~= n)
                        try
                            t4 = this.buildVisitor.filenameT4(stagedImgs{n}, stagedImgs{m});
                            
                            %% CONTINUALLY UPDATE T4

                            this.buildVisitor.align_crossModal( ...
                                'dest',       blurredImgs{m}, ...
                                'source',     blurredImgs{n}, ...
                                'destMask',   maskedImgs{m}, ...
                                'sourceMask', 'none', ...
                                't4',         t4, ...
                                't4img_4dfp', false, ...
                                'useMetricGradient', ipr.useMetricGradient, ...
                                'log',        this.imageRegLog);
                            copyfile(t4, ipr.cacheLoc);
                            
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
            this = this.estimateErr(stagedImgs);
            this.deleteTrash;
        end
        
        %%
         
        function         assessT4Failures(this, t4fails)
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
        function this  = estimateErr(this, stagedImgs)
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
            this.appendVMRBLog([datestr(now) '\n']);
            this.appendVMRBLog( ...
                sprintf('VariableMaskT4ResolveBuilder.estimateErr: stagedImgs->%s\n', ...
                cell2str(stagedImgs, 'AsRow', true)));
            this.appendVMRBLog( ...
                sprintf('VariableMaskT4ResolveBuilder.estimateErr: this.t4_resolve_err->%s\n', ...
                mat2str(this.t4_resolve_err)));            
        end
        function         appendVMRBLog(this, s)
            if (~lexist(this.vmRBLog))
                mlbash(sprintf('touch %s', this.vmRBLog));
            end
            mlbash(sprintf('echo \"%s\" >> %s', s, this.vmRBLog));
            fprintf(s);
        end     
        function fqfps = lazyMasksForImages__(this, ipr, stagedImgs, varargin)
            %% LAZYMASKSFORIMAGES__ has not worked.  Results archived in $PPG/jjlee2/HYGLY25/Vall_2018ma16_*.
            %  @param ipr.maskForImages is 'none' or the fileprefix/name of an anatomical image which will be 
            %  thresholded to generate a mask; the anatomical image must have transverse orientation.   
            %  See also:  mlfourdfp.FourdfpVisitor.transverse_t4.
            %  @param this.theImages must be populated with the images to be masked.
            %  @param stagedImgs is cell-arrayed.
            %  @return fqfps, a cell-array of binary masks, identically sized to this.theImages.
            
            fqfps = cell(1, length(ipr.source));            
            if (~iscell(ipr.maskForImages) && strcmpi(ipr.maskForImages, 'none'))
                fqfps = cellfun(@(x) 'none', fqfps, 'UniformOutput', false);
                return
            end   
            
            % use mask for first source for all other sources
            if (~iscell(ipr.maskForImages))
                ipr.maskForImages = repmat(ensureCell(ipr.maskForImages), size(ipr.source));
            end
            assertSizeEqual(ipr.maskForImages, ipr.source);            
            ii = 1;  
            if (strcmp(ipr.maskForImages{ii}, 'Msktgen'))
                try
                    mg = mlpet.Msktgen('sessionData', this.sessionData);
                    mskt = mg.constructMskt( ...
                        'source', ipr.source{ii}, ...
                        'intermediaryForMask', this.sessionData.T1001, ...
                        'sourceOfMask', fullfile(this.sessionData.sessionPath, 'brainmask.4dfp.hdr'), ...
                        'blurForMask', 33);
                    fqfps{ii} = mskt.fqfileprefix;
                catch ME
                    handexcept(ME);
                end
            end
            for ii = 2:length(ipr.source)
                t4 = this.buildVisitor.filenameT4(stagedImgs{ii}, stagedImgs{1}); % from iteration without mask
                assert(lexist(t4));
                fqfps{ii} = [ipr.source{ii} '_mskt'];
                mg.t4img_4dfp(t4, mskt.fqfileprefix, 'out', fqfps{ii}, 'options', ['-O' ipr.source{1}]);
            end
        end        
        function         purgeT4s(this, ipr)
            stagedImgs = this.lazyStageImages(ipr);
            len = length(stagedImgs);
            for m = 1:len
                for n = 1:len
                    if (m ~= n)
                        try
                            deleteExisting(this.buildVisitor.filenameT4(stagedImgs{n}, stagedImgs{m}));
                        catch ME
                            handexcept(ME);
                        end
                    end
                end
            end
        end
        
 		function this = VariableMaskT4ResolveBuilder(varargin)
 			%% VARIABLEMASKT4RESOLVEBUILDER
 			%  @param .

 			this = this@mlfourdfp.CompositeT4ResolveBuilder(varargin{:});
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

