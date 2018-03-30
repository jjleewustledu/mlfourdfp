classdef (Abstract) AbstractUmapResolveBuilder0 < mlfourdfp.MMRResolveBuilder0
	%% ABSTRACTUMAPRESOLVEBUILDER  

	%  $Revision$
 	%  was created 14-Nov-2016 22:37:22
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.
 	
    
    properties (Constant)
        REPLACE_COMPLETED = false
    end

	properties
        mprFromFdgNacResolveSequence = 3    
        resolveSequenceNRevisions = 2
    end
    
	methods         
 		function this  = AbstractUmapResolveBuilder0(varargin)
 			%% ABSTRACTUMAPRESOLVEBUILDER
 			%  Usage:  this = AbstractUmapResolveBuilder()

 			this = this@mlfourdfp.MMRResolveBuilder0(varargin{:});
            this.NRevisions = 1;
            this.blurArg_ = 1.5; % per Avi, 2016oct25
        end
         
        function resolved        = alignO15NACs(this, o15Nacs)
            fdg  = this.sumTimes( ...
                   this.sessionData.fdgNACResolved('typ', 'fqfp'));
            mpr  = this.mprForAlignment; % on fdg
            umap = this.sessionData.umapSynth('typ', 'fqfp'); % on fdg
            if (isempty(mpr))
                resolved = this.resolveSequence(fdg, o15Nacs{:}, umap);
            else
                resolved = this.resolveSequence(fdg, o15Nacs{:}, umap, mpr);
            end
        end 
        function [this,umaps]    = alignUmapToNACFrames(this, umap)
            
            umapsForNACSelect = this.invertedOrdinateUmaps;
            this.t4sImg( ...
                this.invertedOrdinateT4s, ...
                umap, ...
                umapsForNACSelect, ...
                this.tracerNACRevision_fp_r(this.NRevisions)); % cell, char, cell, char
            
            umaps  = cell(1, length(this.frames));
            ordFr  = 0;
            for fr = 1:length(this.frames)
                
                % assembled selected ordinate frames and unselected frames
                if (this.frames(fr))
                    ordFr = ordFr + 1;
                    this.buildVisitor.move_4dfp( ...
                        umapsForNACSelect{ordFr}, ...
                        this.umapsForNACFrame(fr));
                else
                    % set unselected frames <- this.frame0
                    
                    % tmp = tmpFileprefix;
                    % this.buildVisitor.cropfrac_4dfp( ...
                    %  this.firstCrop, 'umapOnCT_on_fdgv1r1_resolved_sumt', tmp); % this.sessionData.fdgUmapLM('typ', 'fqfp')
                    this.buildVisitor.copy_4dfp(umapsForNACSelect{1}, this.umapsForNACFrame(fr));
                end
                umaps{fr} = this.umapsForNACFrame(fr);
            end
        end 
        function [this,umapsOut] = alignUmapToNACs(this, nacs)
            assert(length(nacs) == length(this.sessionDataCache));
            
            N     = length(nacs);
            Numap = N + 2; % cf. this.alignO15NACs for ordering args to resolveSequence;
                           % typically:  {fdg o15Nacs{:} umap mpr}
                              
            umap     = this.sessionData.umapSynth('typ', 'fqfp');
            umapsOut = cell(1, N);
            
            for o = 1:N
                t4 = this.mapcycle(Numap, o+1);
                [pth,fp] = fileparts(nacs{o});
                [~,umapfp] = fileparts(umap);
                umapsOut{o} = fullfile(pth, sprintf('%s_on_%s', umapfp, fp));
                this.buildVisitor.t4img_4dfp(t4, umap, 'out', umapsOut{o}, 'options', ['-O' umap]);
            end
        end  
        
        function dest  = buildTracerNAC(this, varargin)
            %% BUILDTRACERNAC builds 4dfp formatted NAC images.
            %  See also:  mlfourdfp.FourdfpVisitor.sif_4dfp.
            
            ip = inputParser;
            addOptional(ip, 'sessionData', this.sessionData, @(x) isa(x, 'mlpipeline.ISessionData'));
            parse(ip, varargin{:});            
            sessd = ip.Results.sessionData;
            sessd.attenuationCorrected = false;
            lm      = sessd.tracerLM(      'typ', 'fqfp');
            dest    = sessd.tracerRevision('typ', 'fqfp');
            destLoc = sessd.tracerRevision('typ', 'path');
            
            if (this.buildVisitor.lexist_4dfp(dest))
                return
            end
            if (~this.buildVisitor.lexist_4dfp(lm))
                fprintf('mlfourdfp.AbstractUmapResolveBuilder.buildTracerNAC.buildVisitor.sif_4dfp is building %s\n', lm);
                pwd0 = pwd;
                cd(fileparts(lm));
                this.buildVisitor.sif_4dfp(lm);
                cd(pwd0);
            end
            if (~isdir(destLoc))
                mkdir(destLoc);
            end
            ipr = struct('dest', dest, 'source', lm, 'firstCrop', 0.5, 'rnumber', 0);
            this.cropOrCopy(ipr);
        end
        function this  = convertUmapTo4dfp(this)
            pwd0 = pwd;
            if (~lexist(this.sessionData.tracerUmapLM('typ', '4dfp.img')))
                cd(                             this.sessionData.tracerUmapLM('typ', 'path'));
                this.buildVisitor.IFhdr_to_4dfp(this.sessionData.tracerUmapLM('typ', 'v.hdr'));
                cd(pwd0);
            end
        end
        function this  = convertUmapToE7Format(this, varargin)
            sessd = this.sessionData;
            sessd.rnumber = 1;
            
            ip = inputParser;
            addOptional(ip, 'umap', ...
                fullfile(sessd.tracerNACLocation, ...
                    sprintf('%s_on_%s', sessd.umapSynth('typ', 'fp'), sessd.tracerNACRevision('typ', 'fp'))), ...
                @mlfourdfp.FourdfpVisitor.lexist_4dfp);
            parse(ip, varargin{:});
            umap = ip.Results.umap;
            
            flipped = this.buildVisitor.flip_4dfp('z', umap);
            ic = mlfourd.ImagingContext([flipped '.4dfp.ifh']);
            ic = ic.zoomed(this.zoomForInverseCrop);
            ic.noclobber = false;
            ic.saveas([flipped '.4dfp.ifh']);
            movefile( ...
                sprintf('%s.4dfp.img', flipped), ...
                sprintf('%s.v',        umap), 'f');
            delete(sprintf('%s.4dfp.*', flipped));
            delete(sprintf('%sfz.4dfp.*', umap));
        end
        function this  = convertUmapsToE7Format(this, umaps)
            assert(iscell(umaps));
            for fr = 1:length(umaps)
                this = this.convertUmapToE7Format(umaps{fr});
            end
        end  
        function fqt4  = mapcycle(this, frame0, frameF, varargin)
            switch (this.resolveSequenceNRevisions)
                case 1
                    fqt4 = this.map1cycle(frame0, frameF, varargin{:});
                case 2
                    fqt4 = this.map2cycle(frame0, frameF, varargin{:});
                otherwise
                    error('mlfourdfp:unsupportedSwitchCase', ...
                          'O15UmapResolveBuilder.alignUmapToNACs.resolveSequenceNRevisions -> %i', ...
                          this.resolveSequenceNRevisions);
            end
            assert(lexist(fqt4, 'file'));
        end
        function fqt4  = map1cycle(this, frame0, frameF, varargin)
            ip = inputParser;
            addRequired(ip, 'frame0', @isnumeric);
            addRequired(ip, 'frameF', @isnumeric);
            addOptional(ip, 't4output', sprintf('map1cycle_to_%sr1_frame%i_t4', this.resolveSequenceTag, frameF), @ischar);
            parse(ip, frame0, frameF, varargin{:});
            
            pwd0 = pushd(this.sessionData.tracerNACT4Location);
            fqt4 = this.buildVisitor.t4_mul( ...
                     this.T(frame0, 1), ...
                     this.buildVisitor.t4_inv(this.T(frameF, 1)), basename(ip.Results.t4output));
            fqt4 = fullfile(this.sessionData.tracerNACT4Location, fqt4);
            popd(pwd0);
        end
        function fqt4  = map2cycle(this, frame0, frameF, varargin)
            ip = inputParser;
            addRequired(ip, 'frame0', @isnumeric);
            addRequired(ip, 'frameF', @isnumeric);
            addOptional(ip, 't4output', sprintf('map2cycle_to_%sr1_frame%i_t4', this.resolveSequenceTag, frameF), @ischar);
            parse(ip, frame0, frameF, varargin{:});
            
            b = this.buildVisitor;
            
            % 0 := frame0, F := frameF, r := resolved
            pwd0 = pushd(this.sessionData.tracerNACT4Location);
            T_0r_1 = this.T(frame0, 1);
            T_0r_2 = this.T(frame0, 2);
            T_rF_2 = b.t4_inv(this.T(frameF, 2));
            T_rF_1 = b.t4_inv(this.T(frameF, 1));
            fqt4   = b.t4_mul(T_0r_1, b.t4_mul(T_0r_2, b.t4_mul(T_rF_2, T_rF_1)), basename(ip.Results.t4output));  
            fqt4   = fullfile(this.sessionData.tracerNACT4Location, fqt4);
            popd(pwd0);
        end
        function [this,mprToAtlT4] = ...
                         mpr2atl_4dfp(this)
            mprToAtlT4 = [this.sessionData.mpr('typ', 'fp') '_to_' this.atlas('typ', 'fp') '_t4'];
            if (lexist(mprToAtlT4, 'file'))
                return
            end
            
            this.buildVisitor.mpr2atl_4dfp(this.sessionData.mpr('typ', 'fqfp'), 'options', ['-T' this.atlas('typ', 'fqfp')]);
        end
        function mpr   = mprForAlignment(this)
            f18urb = this.f18UmapResolveBuilder_;
            rsr = f18urb.resolveSequenceResolved('typ', 'fqfp');
            this.buildVisitor.extract_frame_4dfp(rsr, this.mprFromFdgNacResolveSequence);
            mpr = sprintf('%s_on_%s', ...
                  f18urb.sessionData.mpr('typ', 'fp'), ...
                  f18urb.sessionData.fdgNACResolved('typ', 'fp'));
            this.buildVisitor.move_4dfp(sprintf('%s_frame%i', rsr, this.mprFromFdgNacResolveSequence), mpr);
        end 
        function         pushNACUmap(this, varargin)
            ip = inputParser;
            addParameter(ip, 'visits', {'V1' 'V2'}, @iscell);
            parse(ip, varargin{:});
            
            sessd = this.sessionData;
            cd(fullfile(sessd.sessionPath));            
            for iv = 1:length(ip.Results.visits)
                cd(fullfile(sessd.sessionPath, ip.Results.visits{iv}, ...
                    sprintf('%s_%s-NAC', upper(sessd.tracer), ip.Results.visits{iv}), ''));
                try
                    mlbash(sprintf('scp -qr %s %s:%s', ...
                        fullfile( ...
                            sessd.sessionPath, ip.Results.visits{iv}, ...
                            sprintf('%s_%s-NAC/%s_%s-LM-00-umap.4dfp.*', ...
                                upper(sessd.tracer), ip.Results.visits{iv}, upper(sessd.tracer), ip.Results.visits{iv})), ...
                        this.CLUSTER_HOSTNAME, ...
                        fullfile( ...
                            this.CLUSTER_SUBJECTS_DIR, sessd.sessionFolder, ip.Results.visits{iv}, '') )); 
                catch ME
                    handwarning(ME);
                end
            end
        end
        function         reconvertUmapsToE7Format(this)
            pwd0 = pushd(this.sessionData.fdgNACLocation);
            for fr = 1:length(this.frames)
                this.buildVisitor.extract_frame_4dfp(this.umapsForNAC, fr);
                this.convertUmapToE7Format(this.umapsForNACFrame(fr));
                delete([this.umapsForNACFrame(fr) '.4dfp.*']);
                delete([this.umapsForNACFrame(fr) '_flipz.log']);
            end
            popd(pwd0);
        end
        function loc   = resolveSequenceLocation(this, varargin)
            %  @param named tracer is a string identifier.
            %  @param named snumber is the scan number; is numeric.
            %  @param named typ is string identifier:  folder path, fn, fqfn, ...  
            %  See also:  imagingType.
            %  @param named frame is numeric.
            %  @param named rnumber is the revision number; is numeric.
            %  @returns ipr, the struct ip.Results obtained by parse.            
            %  @returns schr, the s-number as a string.
            
            ipr = this.sessionData.iprLocation(varargin{:});
            tag = this.resolveSequenceTag;
            loc = locationType(ipr.typ, ...
                    fullfile(this.sessionData.tracerNACLocation( ...
                        'tracer', this.sessionData.tracer, ...
                        'snumber', this.sessionData.snumber), ...
                    sprintf('%s%s', upper(tag(1)), tag(2:end)), '')); 
        end
        function obj   = resolveSequenceResolved(this, varargin)
            %  @param named tracer is a string identifier.
            %  @param named snumber is the scan number; is numeric.
            %  @param named typ is string identifier:  folder path, fn, fqfn, ...  
            %  See also:  imagingType.
            %  @param named frame is numeric.
            %  @param named rnumber is the revision number; is numeric.
            %  @returns ipr, the struct ip.Results obtained by parse.            
            %  @returns schr, the s-number as a string.
            
            fqfn = fullfile( ...
                this.resolveSequenceLocation, ...
                sprintf('%sr%i_resolved.4dfp.ifh', this.resolveSequenceTag, this.resolveSequenceNRevisions));
            obj  = this.sessionData.fqfileprefixObject(fqfn, varargin{:});
        end
        function fp    = resolveSequenceTag(this)
            tr = lower(this.sessionData.tracer);
            if (strcmp(tr, 'ho'))
                tr = 'oo';
            end
            fp = sprintf('%sResolveSequencev%i', tr, this.sessionData.vnumber);
        end
        function t4    = T(this, varargin)
            ip = inputParser;
            addRequired(ip, 'f1', @isnumeric);
            addRequired(ip, 'r1', @isnumeric);
            addOptional(ip, 'f2', nan, @isnumeric);
            addOptional(ip, 'r2', nan, @isnumeric);
            addParameter(ip, 'fq', false, @islogical);
            parse(ip, varargin{:});
            
            if (isnan(ip.Results.f2))
                t4__ = sprintf('%sr%i_frame%i_to_resolved_t4', ...
                    this.resolveSequenceTag, ip.Results.r1, ip.Results.f1);
            else
                t4__ = sprintf('%sr%i_frame%i_to_%sr%i_frame%i_t4', ...
                    this.resolveSequenceTag, ip.Results.r1, ip.Results.f1, this.resolveSequenceTag, ip.Results.r2, ip.Results.f2);
            end    
            if (~ip.Results.fq)
                t4 = t4__;
            else
                t4 = fullfile(this.sessionData.tracerNACT4Location, t4__);                
            end
        end
        function umaps = umapsForNAC(this)
            umaps = sprintf('%sUmapsForNAC', lower(this.sessionData.tracer));
        end
        function fp    = umapsForNACFrame(this, fr)
            fp = this.frameFileprefix(this.umapsForNAC, fr);
        end
        function [s,r] = viewUmaps(this)
            [s,r] = mlbash(sprintf( ...
                'fslview %s.4dfp.hdr %s.4dfp.hdr', ...
                this.sessionData.tracerNACRevision, ...
                this.umapsForNAC));
        end
        function [s,r] = viewTracerResolved(this)
            [s,r] = mlbash(sprintf( ...
                'fslview %s_resolved_b55.4dfp.hdr', ...
                this.sessionData.tracerNACRevision('typ', 'fqfp')));
        end
    end 
    
    %% PROTECTED
    
    properties (Access = protected)
        carneyUmapBuilder_
        f18UmapResolveBuilder_
    end
    
    methods (Access = protected)
        function resolvedToFramesNACT4s ... 
                          = invertedOrdinateT4s(this)
            t4ss = cell(1, this.NRevisions);
            for r = this.NRevisions:-1:1
                t4ss{r} = this.t4sInv(this.framesToResolvedT4s(this.tracerNACRevision_fp_r(r), this.fortranFrameIndices));
            end
            
            resolvedToFramesNACT4s = t4ss{this.NRevisions};
            for r = this.NRevisions:-1:2
                resolvedToFramesNACT4s = this.t4sMul(resolvedToFramesNACT4s, t4ss{r-1});
            end            
        end
        function umaps    = invertedOrdinateUmaps(this)
            umaps = fullfile(this.sessionData.tracerNACResolved('typ', 'path'), [this.umapsForNAC 'Select']);
            umaps = cellfun( ...
                @(x) this.frameFileprefix(umaps, x), ...
                num2cell(this.fortranFrameIndices), ...
                'UniformOutput', false);            
        end
        function resolved = resolveSequence(this, varargin)
            N = length(varargin);            
            pasted = fullfile(this.sessionData.vLocation, this.resolveSequenceTag);
            this.pasteImages(varargin, pasted);
            this.resolve( ...
                'dest', [pasted '_resolved'], ...
                'source', pasted, ...
                'sourceBlur', this.blurArg * ones(1,N), ...
                'firstCrop', 1, ...
                'frames', ones(1,N), ...
                'NRevisions', this.resolveSequenceNRevisions);
            resolved = sprintf('%sr%i_resolved', pasted, this.resolveSequenceNRevisions);
        end
        function tr       = tracerForNAC(this)
            tr = sprintf('%sForNAC', lower(this.sessionData.tracer));
        end 
        function fp       = tracerNACRevision_fp_r(this, r)
            fp = sprintf('%sv%ir%i', lower(this.sessionData.tracer), this.sessionData.vnumber, r);
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

