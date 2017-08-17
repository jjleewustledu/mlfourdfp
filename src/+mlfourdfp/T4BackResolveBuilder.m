classdef T4BackResolveBuilder < mlfourdfp.AbstractT4ResolveBuilder
	%% T4BACKRESOLVEBUILDER  

	%  $Revision$
 	%  was created 07-Aug-2017 05:41:12 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.2.0.538062 (R2017a) for MACI64.  Copyright 2017 John Joowon Lee.

    
	methods		  
 		function this = T4BackResolveBuilder(varargin)
 			%% T4BACKRESOLVEBUILDER

 			this = this@mlfourdfp.AbstractT4ResolveBuilder(varargin{:});
            ip = inputParser;
            ip.KeepUnmatched = true;           
            addOptional( ip, 'cctor', []);
            addParameter(ip, 'blurArg', 5.5, @isnumeric);
            addParameter(ip, 'indicesLogical', true, @islogical);
            addParameter(ip, 'theImages', {}, @(x) iscell(x) || ischar(x));
            addParameter(ip, 'indexOfReference', 1, @isnumeric);
            parse(ip, varargin{:});
            
            if (isempty(ip.Results.cctor))
                import mlfourdfp.*;
                this.imageComposite_ = ImageFrames(this, ...
                    'indicesLogical', ip.Results.indicesLogical, ...
                    'theImages', FourdfpVisitor.ensureSafeFileprefix(ip.Results.theImages), ...
                    'indexOfReference', ip.Results.indexOfReference);
                this.blurArg_ = ip.Results.blurArg;
            end            
            this.finished = mlpipeline.Finished(this, ...
                'path', this.logPath, ...
                'tag', sprintf('%s_NRev%i_idxOfRef%i', ...
                       lower(this.sessionData.tracerRevision('typ','fp')), this.NRevisions, this.indexOfReference));
            cd(this.sessionData.tracerLocation);
 		end
                
        function this         = resolve(this, varargin)
            %% BACKRESOLVE iteratively calls t4_resolve and writes a log.
            %  @param dest       is a f.q. fileprefix.
            %  @param destMask   "
            %  @param source     "
            %  @param sourceMask "
            %  @param destBlur   is the fwhm blur applied by imgblur_4dfp to dest.            
            %  @param sourceBlur is the fwhm blur applied by imgblur_4dfp to source.
            %  @param maskForImages is a f.q. fileprefix.
            %  @param indicesLogical is logical.
            %  @param t40        is the initial t4-file for the transformation:  transverse is default.
            %  @param resolveTag is char.
            %  @param log        is the f.q. filename of the log file.
            
            import mlfourdfp.*;
            ip = inputParser;
            addParameter(ip, 'dest',           '',                  @ischar);
            addParameter(ip, 'source',         this.theImages,      @FourdfpVisitor.lexist_4dfp);
            addParameter(ip, 'destMask',       'none',              @ischar);
            addParameter(ip, 'sourceMask',     'none',              @ischar);
            addParameter(ip, 'destBlur',       this.blurArg,        @isnumeric); % fwhh/mm
            addParameter(ip, 'sourceBlur',     this.blurArg,        @isnumeric); % fwhh/mm 
            addParameter(ip, 'maskForImages',  'maskForImages',     @ischar);
            addParameter(ip, 'indicesLogical', this.indicesLogical, @islogical);
            addParameter(ip, 't40',            this.buildVisitor.transverse_t4, @(x) ischar(x) || iscell(x));
            addParameter(ip, 'resolveTag',     this.resolveTag,     @ischar);
            addParameter(ip, 'log',            '/dev/null',         @ischar);
            parse(ip, varargin{:});
            this.indicesLogical = ip.Results.indicesLogical;
            this.resolveTag = ip.Results.resolveTag;            
            ipr = ip.Results;
            ipr = this.expandBlurs(ipr);
            
            if (isempty(ipr.dest)); ipr.dest = ipr.source; end
            ipr.resolved = ipr.source; % initialize this.revise   
            %if (this.isfinished)  
            %    this = this.alreadyFinalized(ipr);
            %    return
            %end
            while (this.sessionData.rnumber <= this.NRevisions)
                ipr.source = ipr.resolved;
                ipr.dest   = this.fileprefixRevision(ipr.dest, this.sessionData.rnumber);
                [ipr,this] = this.revise(ipr);
                assert(this.sessionData.rnumber < 10);
            end
            this = this.finalize(ipr);
        end
        function [ipr,this]   = revise(this, ipr)
            this.copySourceToDest(ipr); % crop/copy ipr.source to ipr.dest
            this.resolveLog = loggerFilename( ...
                ipr.dest, 'func', 'T4ResolveBuilder_resolveAndPaste', 'path', this.logPath);
            
            ipr = this.resolveAndPaste(ipr); 
            this.teardownRevision;
            this.sessionData_.rnumber = this.sessionData_.rnumber + 1;
        end
        function [ipr,imgFps] = resolveAndPaste(this, ipr)
            %% RESOLVEANDPASTE - preassign ipr.dest, this.resolveTag, this.indexOfReference as needed.
            %  @param ipr is a struct w/ field dest, a string fileprefix || is a string.
            
            dest_ = mybasename(ipr.dest);
            imgFps = this.fileprefixOfReference(ipr, this.indexOfReference); % initial on ipr.dest
            for f = 1:length(this.indicesLogical)
                if (this.indicesLogical(f) && f ~= this.indexOfReference)
                    %                    fileprefix of frame != this.indexOfReference
                    imgFps = [imgFps ' ' this.fileprefixIndexed(dest_, f)]; %#ok<AGROW>
                end
            end
            
            this.buildVisitor.t4_resolve( ...
                this.resolveTag, imgFps, ...
                'options', '-v -m -s', 'log', this.resolveLog);          
            ipr.resolved = sprintf('%s_%s', dest_, this.resolveTag); 
        end        
        function this         = finalize(this, ipr)            
            this.ipResults_ = ipr;
            this.sessionData.rnumber = this.NRevisions;
            this.product_ = mlpet.PETImagingContext([ipr.resolved '.4dfp.ifh']);              
            this.buildVisitor.imgblur_4dfp(ipr.resolved, this.blurArg);
            this.teardownResolve(ipr);
            this.finished.touchFinishedMarker;          
        end
        function this         = alreadyFinalized(this, ipr)            
            dest = this.fileprefixRevision(ipr.dest, this.NRevisions);
            ipr.resolved = sprintf('%s_%s', dest, this.resolveTag);            
            this.ipResults_ = ipr;
            this.sessionData.rnumber = this.NRevisions;
            this.product_ = mlpet.PETImagingContext([ipr.resolved '.4dfp.ifh']);
        end
        function                teardownResolve(this, ipr)
            if (this.keepForensics); return; end
            
            for r = 1:this.NRevisions                
                fp0 = this.fileprefixRevision(ipr.dest, r);
                for f = 1:length(this.indicesLogical)
                    if (this.indicesLogical(f))
                        %delete(sprintf('%s_frame%i.4dfp.*', fp0, f));
                        delete(sprintf('%s_frame%i_b*.4dfp.*', fp0, f));
                        delete(sprintf('%s_frame%i_C*.4dfp.*', fp0, f));
                        %delete(sprintf('%s_frame%i_%s.4dfp.*', fp0, f, this.resolveTag));
                        delete(sprintf('%s_frame%i_g*.nii.gz', fp0, f));
                    end
                end
                sessd = this.sessionData;
                sessd.rnumber = r;
            end            
            delete(sprintf('%s_*_*.4dfp.*', ipr.maskForImages));
        end
        
        %% UTILITY         
              
        function         copySourceToDest(this, ipr)
            if (1 == this.sessionData.rnumber)
                try
                    if (~this.buildVisitor.lexist_4dfp(ipr.dest) && ...
                        ~strcmp(ipr.source, ipr.dest))
                        this.buildVisitor.copy_4dfp(ipr.source, ipr.dest);
                    end
                catch ME
                    handexcept(ME);
                end
                return
            end
            this.copyResolvedToNewRevision(ipr);
        end
        function         copyResolvedToNewRevision(this, ipr)
            %% COPYRESOLVEDTONEWREVISION opportunistically reuses existing files from the last iteration
            
            this.buildVisitor.copy_4dfp( ...
                this.fileprefixResolved(ipr.dest, this.sessionData.rnumber-1), ...
                this.fileprefixRevision(ipr.dest, this.sessionData.rnumber));
            dt = mlsystem.DirTool( ...
                sprintf('%s_frame*_%s.4dfp.ifh', this.sessionData.tracerRevision('typ','fp'), this.resolveTag));
            for f = 1:length(dt.fns)
                resolvedFrame = mybasename(dt.fns{f});
                idx_   = regexp(resolvedFrame, '_frame\d+');
                idx__  = regexp(resolvedFrame, sprintf('_%s$', this.resolveTag));
                newRev = [ipr.dest resolvedFrame(idx_:idx__-1)];
                try
                    this.buildVisitor.move_4dfp(resolvedFrame, newRev);
                catch ME
                    handexcept(ME);
                end
            end
        end
        function fp    = fileprefixIndexed(~, fp, fr)
            assert(ischar(fp));
            assert(isnumeric(fr));
            fp = sprintf('%s_frame%i', fp, fr);
        end
        function fp    = fileprefixIndexedResolved(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fp', @ischar);
            addRequired(ip, 'fr', @isnumeric);
            addOptional(ip, 'tag', this.resolveTag, @ischar);
            parse(ip, varargin{:});
            fp = sprintf('%s_%s', this.fileprefixIndexed(ip.Results.fp, ip.Results.fr), ip.Results.tag);
        end
        function fp    = fileprefixOfReference(this, ipr, varargin)
            ip = inputParser;
            addOptional(ip, 'indexOfRef', this.indexOfReference, @isnumeric)
            parse(ip, varargin{:});
            
            fp = this.fileprefixIndexed(mybasename(ipr.dest), ip.Results.indexOfRef);
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

