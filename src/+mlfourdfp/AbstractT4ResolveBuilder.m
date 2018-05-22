classdef (Abstract) AbstractT4ResolveBuilder < mlpipeline.AbstractSessionBuilder & mlfourdfp.IT4ResolveBuilder
	%% ABSTRACTT4RESOLVEBUILDER  

	%  $Revision$
 	%  was created 10-Jan-2017 22:33:57
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.  Copyright 2017 John Joowon Lee.
 	

	properties
        imageRegLog
        maskForImagesThreshFactor = 0.25
        NRevisions        
        resolveLog
        skipT4imgAll = false
        t4_resolve_err
    end
    
    properties (Dependent)
        blurArg
        epoch
        epochLabel      
        gaussArg
        imageComposite
        indicesLogical
        indexOfReference
        petBlur
        referenceImage
        referenceWeight
        resolveTag
        sourceImage
        sourceWeight
        t4s
        theImages
    end
    
    methods (Static)
        function       diaryv(name)
            assert(ischar(name));
            if (~isempty(getenv('PRINTV')))
                diary(sprintf('mlfourdfp_AbstractT4ResolveBuilder%s.txt', name)); 
            end
        end
        function en  = ensureLocalFourdfp(en)
            en_ = en;
            en  = mybasename(mlfourdfp.FourdfpVisitor.ensureLocalFourdfp(en));
            if (iscell(en_))
                en = ensureCell(en);
            end
        end
        function obj = ensureSumtSaved(obj, varargin)
            %  @param named typ; see also mlfourd.ImagingContext.imagingType.
            %  @return fqfp of sumt after saving (default).
            
            ip = inputParser;
            addParameter(ip, 'typ', 'fqfp', @ischar);
            parse(ip, varargin{:});
            
            ic = mlfourd.ImagingContext(obj);
            ic = ic.timeSummed;
            ffp = ic.fourdfp;
            ffp.save;
            obj = ic.imagingType(ip.Results.typ, ic);
        end
        function f   = frameNumber(str, offset)
            names = regexp(str, '\w+(-|_)(F|f)rame(?<f>\d+)', 'names');
            f = str2double(names.f) + offset;
        end
        function       printv(varargin)
            if (~isempty(getenv('PRINTV')))
                fprintf('mlfourdfp.AbstractT4ResolveBuilder.');
                fprintf(varargin{:});
                fprintf('\n');
            end
        end
    end
    
    methods 
        
        %% GET/SET
        
        function g    = get.blurArg(this)
            if (~isempty(this.blurArg_))
                g = this.blurArg_;
                return
            end
            
            %g = 5.5; % Avi's empirical suggestion for MR, PET alignment
            g = this.petBlur * 1.2790697674; 
        end
        function this = set.blurArg(this, s)
            assert(isnumeric(s));
            this.blurArg_ = s;
        end
        function g    = get.epoch(this)
            g = this.sessionData_.epoch;
        end
        function this = set.epoch(this, s)
            assert(isnumeric(s));
            this.sessionData_.epoch = s;
        end
        function g    = get.epochLabel(this)
            g = this.sessionData_.epochLabel;
        end
        function g    = get.gaussArg(this)
            %g = 1.1;
            g = this.F_HALF_x_FWHM/this.blurArg;
        end
        function g    = get.imageComposite(this)
            g = this.imageComposite_;
        end
        function this = set.imageComposite(this, s)
            assert(isa(s, 'mlfourdfp.IImageComposite'));
            this.imageComposite_ = s;
        end
        function g    = get.indexOfReference(this)
            g = this.imageComposite.indexOfReference;
        end
        function this = set.indexOfReference(this, s)
            this.imageComposite_.indexOfReference = s;
        end
        function g    = get.indicesLogical(this)
            g = this.imageComposite.indicesLogical;
        end
        function this = set.indicesLogical(this, s)
            this.imageComposite.indicesLogical = s;
        end
        function g    = get.petBlur(this)
            g = this.sessionData.petBlur;
        end
        function g    = get.referenceWeight(~)
            g = [];
        end
        function g    = get.referenceImage(this)
            g = this.imageComposite.referenceImage;
        end
        function g    = get.resolveTag(this)
            g = this.sessionData.resolveTag;
        end
        function this = set.resolveTag(this, s)
            assert(ischar(s));
            this.sessionData_.resolveTag = s;
        end
        function g    = get.sourceImage(this)
            g = this.imageComposite.sourceImage;
        end
        function g    = get.sourceWeight(~)
            g = [];
        end
        function g    = get.t4s(this)
            g = this.t4s_;
        end
        function g    = get.theImages(this)
            if (isempty(this.imageComposite))
                g = [];
                return
            end
            g = this.imageComposite.theImages;
        end
        function this = set.theImages(this, s)
            this.imageComposite_.theImages = s;
        end
    
        %%
        
        function this  = alreadyFinalized(this, ipr)
            dest         = this.fileprefixRevision(ipr.dest, this.NRevisions);
            ipr.resolved = sprintf('%s_%s', dest, this.resolveTag);
            this         = this.buildProduct(ipr);
        end
        function         appendImageRegLog(this, s)
            if (~lexist(this.imageRegLog))
                mlbash(sprintf('touch %s', this.imageRegLog));
            end
            mlbash(sprintf('echo \"%s\" >> %s', s, this.imageRegLog));
            fprintf(s);
        end
        function fp    = buildMaskForImage(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fp1', @(x) lexist([x '.4dfp.ifh']));
            addOptional(ip, 'fp',  ['mask_' varargin{1}], @ischar);
            parse(ip, varargin{:});
            
            import mlfourd.* mlfourdfp.*;
            fqfp1 = [ip.Results.fp1 '_g0_1'];
            bv = this.buildVisitor;
            if (~bv.lexist_4dfp(fqfp1))
                fqfp_ = bv.gauss_4dfp(ip.Results.fp1, 0.1);
                fqfp1 = strrep(fqfp_, '0.1','0_1');
                bv.move_4dfp(fqfp_, fqfp1);
            end
            ic  = ImagingContext([fqfp1 '.4dfp.ifh']);
            ic  = ic.thresh(ic.numericalNiftid.dipmax*this.maskForImagesThreshFactor);
            ic  = ic.binarized;
            ic.saveas([ip.Results.fp '.4dfp.ifh']);
            fp  = ic.fqfileprefix;
        end
        function fp    = buildMaskForImages(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fp1', @(x) lexist([x '.4dfp.ifh']));
            addRequired(ip, 'fp2', @(x) lexist([x '.4dfp.ifh']));
            addOptional(ip, 'fp',  ['mask_' varargin{1} '_' varargin{2}], @ischar);
            parse(ip, varargin{:});
            
            import mlfourd.* mlfourdfp.*;
            fqfp1 = [ip.Results.fp1 '_g0_1'];
            fqfp2 = [ip.Results.fp2 '_g0_1'];
            bv = this.buildVisitor;
            if (~bv.lexist_4dfp(fqfp1))
                fqfp_ = bv.gauss_4dfp(ip.Results.fp1, 0.1);
                fqfp1 = strrep(fqfp_, '0.1','0_1');
                bv.move_4dfp(fqfp_, fqfp1);
            end
            if (~bv.lexist_4dfp(fqfp2))
                fqfp_ = bv.gauss_4dfp(ip.Results.fp2, 0.1);
                fqfp2 = strrep(fqfp_, '0.1','0_1');
                bv.move_4dfp(fqfp_, fqfp2);
            end
            ic1 = ImagingContext([fqfp1 '.4dfp.ifh']);
            ic2 = ImagingContext([fqfp2 '.4dfp.ifh']);
            ic  = ImagingContext(ic1.numericalNiftid + ic2.numericalNiftid);
            ic  = ic.thresh(ic.numericalNiftid.dipmax*this.maskForImagesThreshFactor);
            ic  = ic.binarized;
            ic.saveas([ip.Results.fp '.4dfp.ifh']);
            fp  = ic.fqfileprefix;
        end   
        function ipr   = expandBlurs(this, ipr)
            %  @return ipr.sourceBlur has size(this.indicesLogical)
            %  @return ipr.destBlur   has size(this.indicesLogical)
            
            assert(isstruct(ipr));
            assert(isfield( ipr, 'sourceBlur'));
            assert(isfield( ipr, 'destBlur'));
            ipr.sourceBlur = ipr.sourceBlur .* this.indicesLogical;
            ipr.destBlur   = ipr.destBlur   .* this.indicesLogical;
        end
        function ipr   = expandMasks(~, ipr)
            %  @return ipr.sourceMask is cell with size(this.indicesLogical)
            %  @return ipr.destMask   is cell with size(this.indicesLogical)
            
            assert(isstruct(ipr));
            assert(isfield( ipr, 'sourceMask'));
            assert(isfield( ipr, 'destMask'));
            if (~iscell(ipr.sourceMask))
                assert(ischar(ipr.sourceMask));
                ipr.sourceMask = repmat({ipr.sourceMask}, size(ipr.indicesLogical));
            end
            if (~iscell(ipr.destMask))
                assert(ischar(ipr.destMask));
                ipr.destMask = repmat({ipr.destMask}, size(ipr.indicesLogical));
            end
        end
        function fqfp  = fileprefixBlurred(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fqfp', @ischar);
            addOptional(ip, 'blur', this.blurArg, @isnumeric);
            parse(ip, varargin{:});
            fqfp = this.buildVisitor.fileprefixBlurred(ip.Results.fqfp, ip.Results.blur);
        end        
        function fqfp  = fileprefixGaussed(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fqfp', @ischar);
            addOptional(ip, 'gauss', this.gaussArg, @isnumeric);
            parse(ip, varargin{:});
            fqfp = this.buildVisitor.fileprefixGaussed(ip.Results.fqfp, ip.Results.gauss);
        end
        function fp    = fileprefixIndexed(~, ipr)
            assert(isfield(ipr, 'dest'));
            assert(isfield(ipr, 'currentIndex'));
            fp = sprintf('%s_frame%i', ipr.dest, ipr.currentIndex);
        end
        function fp    = fileprefixResolved(this, fp, rnumber)
            fp = sprintf('%s_%s', this.fileprefixRevision(fp, rnumber), this.resolveTag);
        end   
        function fp    = fileprefixRevision(this, fp, rnumber)
            %% FILEPREFIXREVISION strips this.resolveTag and r[0-9] from fp; replaces r[0-9].
            %  @param fp, e.g., fdgv1r1_resolved
            %  @param rnumber is numeric
            %  @returns fp, e.g., fdgv1r2
            
            assert(ischar(fp));
            assert(isnumeric(rnumber));
            
            idx = regexp(fp, sprintf('_%s$', this.resolveTag));
            if (~isempty(idx))
                fp = fp(1:idx-1);
            end
            if (length(fp) > 2)
                if (~isempty(regexp(fp(end-2:end), 'r[0-9]', 'match')))
                    fp = fp(1:end-2);
                end
            end
            fp = sprintf('%sr%i', fp, rnumber);
        end  
        function fqfp  = fileprefixSumt(this, fqfp)
            assert(ischar(fqfp));
            fqfp = sprintf('%s_sumt%i', fqfp, sum(this.indicesLogical));
        end
        function this  = finalize(this, ipr)
            this.t4ResolveError_.logger.save;
            this = this.buildProduct(ipr);            
            this.teardownResolve(ipr);
            this.finished.touchFinishedMarker;  
        end
        function this  = img2atl(this, fqfp)
            sd = this.sessionData;
            imgToAtlT4 = [fqfp '_to_' sd.atlas('typ', 'fp') '_t4'];
            if (~lexist(imgToAtlT4, 'file'))
                pwd_ = pushd(fileparts(fqfp));
                this.buildVisitor.mpr2atl_4dfp(fqfp, 'options', ['-crossmodal -T' sd.atlas('typ', 'fqfp')]);
                popd(pwd_);
            end
            if (~lexist(mybasename(imgToAtlT4), 'file'))
                this.buildVisitor.lns(imgToAtlT4);
            end
        end
        function fqfp  = lazyBlurImage(this, ipr, blur)
            %% LAZYBLURIMAGE uses specifiers in ipr; will not replace any existing image
            %  @param ipr is a struct
            %  @param blur is numeric
            
            fqfp_ = this.fileprefixIndexed(ipr);
            fqfp  = this.fileprefixBlurred(fqfp_, blur);
            if (~this.buildVisitor.lexist_4dfp(fqfp))
                this.buildVisitor.imgblur_4dfp(fqfp_, blur);
            end
        end
        function fqfps = lazyBlurImages(this, ipr)
            %% LAZYBLURIMAGES uses specifiers in ipr; will not replace any existing indicesLogical
            %  @param ipr is a struct
            
            fqfps = {};
            for il = 1:length(this.indicesLogical)
                if (this.indicesLogical(il))
                    ipr.currentIndex = il;
                    fqfps = [fqfps this.lazyBlurImage(ipr, ipr.destBlur(il))]; %#ok<AGROW>
                end
            end
            assert(length(fqfps) == sum(this.indicesLogical));
        end    
        function msk   = maskBoundaries(this, fqfp)
            this.buildVisitor.nifti_4dfp_ng(fqfp);
            ic  = mlfourd.ImagingContext(fqfp);
            ic  = ic.ones;
            ic.noclobber = false;
            ic.saveas([fqfp '_ones']);
            this.buildVisitor.nifti_4dfp_4([fqfp '_ones']);
            msk = this.zeroSlicesOnBoundaries([fqfp '_ones'], 3);
        end           
        function this  = mpr2atl(this)
            sd = this.sessionData;
            this.mprToAtlT4 = [sd.mpr('typ', 'fqfp') '_to_' sd.atlas('typ', 'fp') '_t4'];
            if (~lexist(this.mprToAtlT4, 'file'))
                pwd_ = pushd(sd.mpr('typ', 'path'));
                this.buildVisitor.mpr2atl_4dfp(sd.mpr('typ', 'fp'), 'options', ['-T' sd.atlas('typ', 'fqfp')]);
                popd(pwd_);
            end
            if (~lexist(mybasename(this.mprToAtlT4), 'file'))
                this.buildVisitor.lns(this.mprToAtlT4);
            end
        end
        function this  = msktgenImg(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fqfp', @lexist_4dfp);
            addOptional(ip, 'atl', fullfile(getenv('REFDIR'), 'TRIO_Y_NDC'), @lexist_4dfp);
            parse(ip, varargin{:});
            fqfp = ip.Results.fqfp;
            atl  = ip.Results.atl;
            
            log = sprintf('msktgenImg_%s.log', datestr(now, 30));
            this.buildVisitor.mpr2atl_4dfp(fqfp, 'options', ['-crossmodal -T' atl], 'log', log);
            this.buildVisitor.msktgen_4dfp(fqfp, 'options', ['-T' atl], 'log', log);
        end 
        function this  = msktgenMprage(this, varargin)
            this = this.buildVisitor.msktgenMprage(this, varargin{:});
        end   
        function pth   = onAtlasPath(this)
            pth = fullfile(this.sessionData.tracerLocation, 'Atlas', '');
            %ensuredir(pth); % CRUFT?
        end
        function r     = resolvePair(this, f1, f2)
            [~,r] = this.buildVisitor.t4_resolve( ...
                this.resolveTag, [f1 ' ' f2], ...
                'options', '-v -m -s');
        end
        function tag   = resolveTagFrame(this, varargin)
            tag = this.sessionData.resolveTagFrame(varargin{:});
        end    
        function fqfp  = sumTimes(this, varargin)
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'fqfpDyn', @lexist_4dfp);
            addParameter(ip, 'indicesLogical', this.indicesLogical, @isnumeric);
            parse(ip, varargin{:});
            
            if (sum(ip.Results.indicesLogical) < this.imageComposite.length)
                fqfp = this.sumTimesImagingContext(varargin{:});
            else
                fqfp = this.sumTimesActmapf4dfp(varargin{:});
            end
        end
        function [rmsdeg,rmsmm] = t4_resolve_errParser(~, r)
            %% T4_RESOLVE_PAIR
            %  @param r is the result from a call to mlbash.
            %  @return rmsdeg is numeric error found by t4_resolve for pair.
            %  @return rmsmm  is numeric error found by t4_resolve for pair.
            % 
            % jjlee@william
            % [ /data/nil-bluearc/raichle/PPGdata/jjlee2/HYGLY14/V1 ]$ t4_resolve T1001 fdgv1r2_op_fdgv1e1to4r1_frame4_sumt -oT1001_to_op_fdgv1r2_op_fdgv1e1to4r1_frame4_sumt
            % $Id: t4_resolve.c,v 1.6 2013/09/12 00:31:12 avi Exp $
            % n=2
            % rms error averaged over   2 pairs
            %     0.2726    0.1902    0.2722    0.0986    0.2615    0.1313
            % pairs total rotation error        0.4296 (rms deg)
            % pairs total translation error     0.3088 (rms mm)
            % t4_resolve: VOI rms radius=50.0000
            % t4_resolve: begin Gauss-Newton trajectory estimation
            % estimated observational error based on   2 observations
            %     0.3722    0.2711    0.3963    0.0046    0.4277    0.2680
            % estimate total rotation error     0.6075 (rms deg)
            % estimate total translation error  0.5047 (rms mm)
            % rigid body rmserr       0.2113
            % JtJ condition number0.1638E+05
            % scale rmserr            0.0000
            % err=0.535842 err0=0.000000
            % estimated observational error based on   2 observations
            %     0.2845    0.1895    0.2715    0.0639    0.2449    0.0563
            % estimate total rotation error     0.4365 (rms deg)
            % estimate total translation error  0.2593 (rms mm)
            % rigid body rmserr       0.1330
            % JtJ condition number0.1633E+05
            % scale rmserr            0.0000
            % err=0.212365 err0=0.535842
            % estimated observational error based on   2 observations
            %     0.2846    0.1895    0.2715    0.0633    0.2449    0.0561
            % estimate total rotation error     0.4366 (rms deg)
            % estimate total translation error  0.2591 (rms mm)
            % rigid body rmserr       0.1330
            % JtJ condition number0.1633E+05
            % scale rmserr            0.0000
            % err=0.212343 err0=0.212365
            % estimated observational error based on   2 observations
            %     0.2846    0.1895    0.2715    0.0633    0.2449    0.0561
            % estimate total rotation error     0.4366 (rms deg)
            % estimate total translation error  0.2591 (rms mm)
            % rigid body rmserr       0.1330
            % JtJ condition number0.1633E+05
            % scale rmserr            0.0000
            % err=0.212343 err0=0.212343
            % estimated observational error based on   2 observations
            %     0.2846    0.1895    0.2715    0.0633    0.2449    0.0561
            % estimate total rotation error     0.4366 (rms deg)
            % estimate total translation error  0.2591 (rms mm)
            % rigid body rmserr       0.1330
            % JtJ condition number0.1633E+05
            % scale rmserr            0.0000
            % err=0.212343 err0=0.212343
            % Writing: T1001_to_op_fdgv1r2_op_fdgv1e1to4r1_frame4_sumt.sub
            % T1001	t4=T1001_to_T1001_to_op_fdgv1r2_op_fdgv1e1to4r1_frame4_sumt_t4
            % fdgv1r2_op_fdgv1e1to4r1_frame4_sumt	t4=fdgv1r2_op_fdgv1e1to4r1_frame4_sumt_to_T1001_to_op_fdgv1r2_op_fdgv1e1to4r1_frame4_sumt_t4

            rmsdeg = nan;
            rmsmm  = nan;
            
            assert(ischar(r) || isstring(r))
            r = splitlines(r);
            %r(r == "") = [];
            r = strip(r);
            foundLabel = strfind(r, 'rms error averaged over');
            ir = 0;
            while (ir <= length(foundLabel))
                ir = ir + 1;
                if (~isempty(foundLabel{ir}))
                    break
                end
            end
            if (ir >= length(foundLabel))
                return
            end
            if (contains(r(ir + 2), 'pairs total rotation error'))
                toknames = regexp(r(ir + 2), '(?<rmsdeg>\d+.\d+) \(rms deg\)', 'names');
                if (~isempty(toknames))
                    rmsdeg = str2double(toknames.rmsdeg);
                end
            end
            if (contains(r(ir + 3), 'pairs total translation error'))
                toknames = regexp(r(ir + 3), '(?<rmsmm>\d+.\d+) \(rms mm\)', 'names');
                if (~isempty(toknames))
                    rmsmm = str2double(toknames.rmsmm);
                end
            end 
        end        
        function err   = t4_resolve_errPairParser(this, obj1, obj2)
            [rmsdeg,rmsmm] = this.t4_resolve_errParser(this.resolvePair(obj1, obj2));
            err = mean([this.rmsarc(rmsdeg) rmsmm], 'omitnan');
        end
        function arc   = rmsarc(~, rmsdeg)            
            arc = 50*pi*rmsdeg/180; % arc at 50 mm from center of mass
        end
        function this  = t4img_4dfp(this, t4fn, source, varargin)
            %% T4IMG_4DFP is a strategy for this.NRevisions == 1 or 2.
            %  @param t4fn is proposed filename structure of the t4 file, e.g., t4source[r2]_to_t4dest_t4.  
            %  @param source is the filename of the new source to be affine transformed according to t4fn;
            %  the present working directory is extracted from source.  
            %  @param named out is the filename of the affine transformed image; filepath := fileparts(source).
            %  @param named ref is the filename of the image that specifies voxel metrics.
            %  @param named options explicitly sends options to mlfourdfp.FourdfpVisitor.t4img_4dfp.  
            %  @returns this with this.product := mlfourd.ImagingContext(out).
            %  N.B.:  t4fn and source may both be cell-arrays.
            
            % recursive
            if (iscell(t4fn) && iscell(source))
                p = cell(size(source));
                for s = 1:length(source)
                    that = this.t4img_4dfp(t4fn{s}, source{s}, varargin{:}, 'out', '');
                    p{s} = that.product_;
                end
                this.product_ = p;
                return
            end
            
            % base-case
            switch (this.NRevisions)
                case 1
                    this = this.t4img_4dfpr1(t4fn, source, varargin{:});
                case 2
                    this = this.t4img_4dfpr2(t4fn, source, varargin{:});
                otherwise
                    error('mlfourdfp:unsupportedSwitchCase', ...
                          'AbstractT4ResolveBuilder.t4img_4dfp.NRevisions->%i', this.NRevisions);
            end
        end
        function this  = t4img_4dfpr1(this, varargin)
            %% T4IMG_4DFPR1 supplies transformations for NRevisions = 1.
            %  @param t4fn is proposed filename structure of the t4 file, e.g., t4source[r2]_to_t4dest_t4.  
            %  @param source is the filename of the new source to be affine transformed according to t4fn;
            %  the present working directory is extracted from source.
            %  @param named out is the filename of the affine transformed image.
            %  @param named ref is the filename of the image that specifies voxel metrics.
            %  @param named options explicitly sends options to mlfourdfp.FourdfpVisitor.t4img_4dfp.  
            %  @returns this with this.product := mlfourd.ImagingContext(out)
            
            ip = inputParser;
            addRequired( ip, 't4fn', @ischar);
            addRequired( ip, 'source', @lexist_4dfp);
            addParameter(ip, 'out', '', @ischar);
            addParameter(ip, 'ref', '', @ischar); % fed to options -Oreference
            addParameter(ip, 'options', '', @ischar); % supplies -Oreference
            parse(ip, varargin{:});  
            
            [t4source,t4dest] = this.buildVisitor.parseFilenameT4(ip.Results.t4fn);  
            t4source = this.markLastRevisionMarking(t4source);
            sourceFqfp = this.markLastRevisionMarking(myfileprefix(ip.Results.source));
            outFqfp = myfileprefix(ip.Results.out);
            if (isempty(outFqfp))
                outFqfp = sprintf('%s_op_%s', this.ensureLastRnumber(sourceFqfp, 1), this.resolveDest(t4dest));
            end
            refFqfp = myfileprefix(ip.Results.ref);
            if (isempty(refFqfp))
                refFqfp = fullfile(fileparts(t4source), this.resolveDest(t4dest));
            end
            options = ip.Results.options;
            if (isempty(options))
                assert(~isempty(refFqfp));
                options = ['-O' refFqfp];
            end
            
            this.buildVisitor.t4img_4dfp( ...
                sprintf('%s_to_%s_t4', this.ensureLastRnumber(t4source,1), t4dest), ...
                this.ensureLastRnumber(sourceFqfp,1), ...
                'out', outFqfp, ...
                'options', options); 
            this.product_ = mlfourd.ImagingContext([outFqfp '.4dfp.ifh']);
        end
        function this  = t4img_4dfpr2(this, varargin)
            %% T4IMG_4DFPR2 supplies transformations for NRevisions = 2.
            %  @param t4fn is proposed filename structure of the t4 file, e.g., t4source[r2]_to_t4dest_t4.  
            %  @param source is the filename of the new source to be affine transformed according to t4fn;
            %  the present working directory is extracted from source.
            %  @param named out is the filename of the affine transformed image; filepath := fileparts(source).
            %  @param named ref is the filename of the image that specifies voxel metrics.
            %  @param named options explicitly sends options to mlfourdfp.FourdfpVisitor.t4img_4dfp.  
            %  @returns this with this.product := mlfourd.ImagingContext(out)
            
            ip = inputParser;
            addRequired( ip, 't4fn', @ischar);
            addRequired( ip, 'source', @lexist_4dfp);
            addParameter(ip, 'out', '', @ischar);
            addParameter(ip, 'ref', '', @ischar); % fed to options -Oreference
            addParameter(ip, 'options', '', @ischar); % supplies -Oreference
            parse(ip, varargin{:});  
            
            [t4source,t4dest] = this.buildVisitor.parseFilenameT4(ip.Results.t4fn);  
                % t4source->${E}/fdgv1${e}r0_frame${e}
                % t4dest->op_fdgv1${e}r1_frame${idxRef}
            t4source = this.markLastRevisionMarking(t4source);
                % t4source->${E}/fdgv1${e}r0_frame${e}
            sourceFqfp = this.markLastRevisionMarking(myfileprefix(ip.Results.source));
                % sourceFqfp ~ E1to9/umapSynth_op_fdgv1e1to9r1_frame${e}
            outFqfp = myfileprefix(ip.Results.out);
                % outFqfp ~ ${E}/umapSynth_op_fdgv1${e}r1_frame${idxRef}
            if (isempty(outFqfp))
                outFqfp = sprintf('%s_op_%s', this.ensureLastRnumber(sourceFqfp, 2), this.resolveDest(t4dest));
            end
            refFqfp = myfileprefix(ip.Results.ref);
                % refFqfp ~ ${E}/fdgv1${e}r1_frame${idxRef}
            if (isempty(refFqfp))
                refFqfp = fullfile(fileparts(t4source), this.resolveDest(t4dest));
            end
            options = ip.Results.options;
                % options ~ -O${E}/fdgv1${e}r1_frame${idxRef}
            if (isempty(options))
                assert(~isempty(refFqfp));
                options = ['-O' refFqfp];
            end                     
            
            assert(lexist(sprintf('%s_to_%s_t4', this.ensureLastRnumber(t4source,1), t4dest), 'file'), ...
                'mlfourdfp:IOErr:fileNotFound', 'AbstractT4ResolveBuilder.t4img_4dfpr2 could not find %s', ...
                sprintf('%s_to_%s_t4', this.ensureLastRnumber(t4source,1), t4dest));          
            assert(lexist(sprintf('%s_to_%s_t4', this.ensureLastRnumber(t4source,2), t4dest), 'file'), ...
                'mlfourdfp:IOErr:fileNotFound', 'AbstractT4ResolveBuilder.t4img_4dfpr2 could not find %s', ...
                sprintf('%s_to_%s_t4', this.ensureLastRnumber(t4source,2), t4dest));
            this.buildVisitor.t4_mul( ...
                sprintf('%s_to_%s_t4', this.ensureLastRnumber(t4source,1),    t4dest), ...
                sprintf('%s_to_%s_t4', this.ensureLastRnumber(t4source,2),    t4dest), ...
                sprintf('%s_to_%s_t4', this.ensureLastRnumbers(t4source,1,2), t4dest));        
            assert(lexist(sprintf('%s_to_%s_t4', this.ensureLastRnumbers(t4source,1,2), t4dest), 'file'), ...
                'mlfourdfp:IOErr:fileNotFound', 'AbstractT4ResolveBuilder.t4img_4dfpr2 could not find %s', ...
                sprintf('%s_to_%s_t4', this.ensureLastRnumbers(t4source,1,2), t4dest));  
            
            this.buildVisitor.t4img_4dfp( ...
                sprintf('%s_to_%s_t4', this.ensureLastRnumbers(t4source,1,2), t4dest), ...
                sourceFqfp, ...
                'out', outFqfp, ...
                'options', options); 
                % t4->${E}/fdgv1${e}r1r2_frame${e}_to_op_fdgv1${e}r1_frame${idxRef}_t4
            fprintf('t4img_4dfpr2 wrote out:\n    %s\n\n', outFqfp);
            this.product_ = mlfourd.ImagingContext([outFqfp '.4dfp.ifh']);
        end
        function pth   = t4Path(this)
            pth = fullfile(this.sessionData.tracerLocation, 'T4', '');
            ensuredir(pth);
        end   
        function         teardownLogs(this)
            ensuredir(this.getLogPath);
            try
                movefiles('*.log', this.getLogPath); 
                movefiles('*.txt', this.getLogPath);   
                movefiles('*.lst', this.getLogPath);    
                movefiles('*.mat0', this.getLogPath);   
                movefiles('*.sub', this.getLogPath); 
            catch ME
                handwarning(ME);
            end
        end
        function         teardownT4s(this)
            if (this.keepForensics); return; end
            
            %% The following will break t4_resolve operations requiring *_t4 files.
            %     try
            %         ensuredir(this.t4Path);
            %         movefiles('*_t4', this.t4Path);
            %         sessd = this.sessionData;
            %         movefile( ...
            %             fullfile(this.t4Path, ...
            %                 sprintf('%s_to_%s_t4', sessd.mpr('typ', 'fp'), sessd.atlas('typ', 'fp'))), ...
            %             pwd);
            %     catch ME
            %         handwarning(ME);
            %     end
        end
        function         teardownRevision(this, ipr)
            if (iscell(ipr.maskForImages))
                return
            end
            try
                deleteExisting(sprintf('%s.4dfp.*', ipr.maskForImages));
                if (this.keepForensics); return; end

                this.teardownLogs;
                this.teardownT4s;
            catch ME
                handwarning(ME);
            end
        end
        function this  = updateFinished(this, varargin)
            %% UPDATEFINISHED, the protected superclass property which is an mlpipeline.Finished
            %  @param tag containing information such as this.sessionData.tracerRevision, this.NRevisions.
            %  this.indexOfReference.
            %  @param tag2.
            %  @param neverTouchFinishfile is boolean.
            %  @param ignoreFinishfile is boolean.
            %  @return property this.finished instantiated with path, tags, the booleans.
            
            this = updateFinished@mlpipeline.AbstractDataBuilder(this, varargin{:});
            this.finished_.tag = sprintf('%s_NRev%i_idxOfRef%i', ...
                this.finished_.tag, this.NRevisions, this.indexOfReference);
        end
        
		function this = AbstractT4ResolveBuilder(varargin)
 			%% ABSTRACTT4RESOLVEBUILDER
            %  @param named NRevisions
            %  @param named keepForensics
            %  @param named resolveTag
            %  @param named theImages
            %  @param named ipResults

 			
            this = this@mlpipeline.AbstractSessionBuilder(varargin{:});
            if (0 == nargin); return; end
            
            %% invoke copy-ctor
            
            if (1 == nargin && isa(varargin{1}, 'mlfourdfp.AbstractT4ResolveBuilder'))
                aCopy = varargin{1};
                this.blurArg_ = aCopy.blurArg_;
                this.imageRegLog = aCopy.imageRegLog;
                this.mprToAtlT4 = aCopy.mprToAtlT4;
                this.msktgenThresh = aCopy.msktgenThresh;
                this.NRevisions = aCopy.NRevisions;
                this.resolveLog = aCopy.resolveLog;
                this.imageComposite_ = aCopy.imageComposite_;
                this.trash_ = aCopy.trash_;
                this.xfm_ = aCopy.xfm_; 
                this.ipResults_ = aCopy.ipResults_;
                return
            end
            
            %% manage parameters 
            
            import mlfourdfp.*;
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'ipResults',     struct([]), @isstruct);
            addParameter(ip, 'keepForensics', false,      @islogical);
            addParameter(ip, 'maskForImages', 'Msktgen',  @(x) ~isempty(x));
            addParameter(ip, 'NRevisions',    2,          @isnumeric);
            addParameter(ip, 'resolveTag',    this.sessionData.resolveTag, @ischar);
            addParameter(ip, 'theImages',     '',         @(x) iscell(x) || ischar(x));
            parse(ip, varargin{:});            
            this.ipResults_     = ip.Results.ipResults;
            this.keepForensics  = ip.Results.keepForensics; % override mlpipeline.AbstractDataBuilder
            this.maskForImages_ = ip.Results.maskForImages;
            this.NRevisions     = ip.Results.NRevisions;
            this.resolveTag     = ip.Results.resolveTag;
            this.theImages      = this.ensureSafeFileprefix(ip.Results.theImages);
            
            %this.logger_.filepath = fullfile(this.sessionData.tracerLocation, 'Log', '');
            %ensuredir(this.logger_.filepath);
            
            this.t4ResolveError_ = T4ResolveError( ...
                'sessionData', this.sessionData, 'logPath', this.getLogPath);
        end     
 	end 

    %% PROTECTED
    
    properties (Access = protected)
        blurArg_      
        imageComposite_
        maskForImages_
        mprToAtlT4
        msktgenThresh = 0
        trash_ = {};
        t4ResolveError_
        t4s_
        xfm_
    end
    
    methods (Access = protected)
        function this = buildProduct(this, ipr)
            this.ipResults_ = ipr;
            this.rnumber = this.NRevisions;            
            this.product_ = mlfourd.ImagingContext([ipr.resolved '.4dfp.ifh']);
        end
        function this = cacheT4s(this, imgFpsc)
            %  @return this.t4s_{r}{1} for r-number r is the reference; size(this.t4s_) == this.NRevisions;
            %  size(this.t4s_{r}) == size(this.indicesLogical).
            
            this.t4s_{this.rnumber} = cell(size(imgFpsc));
            for f = 1:length(imgFpsc)
                try
                    this.t4s_{this.rnumber}{f} = sprintf('%s_to_%s_t4', imgFpsc{f}, this.resolveTag);
                catch ME
                    dispwarning(ME);
                end
            end
        end
        function fp   = clipLastRevisionMarking(~, fp)
            pos = regexp(fp, 'r\d$', 'ONCE');
            if (~isempty(pos))
                fp = fp(1:pos-1);
            end
        end
        function this = copySourceToResolved(this, ipr)
            ipr.dest = this.fileprefixRevision(ipr.dest, this.NRevisions);            
            ipr.resolved = sprintf('%s_%s', ipr.dest, this.resolveTag);
            this.buildVisitor_.copyfile_4dfp(ipr.source, ipr.resolved);
        end
        function this = deleteTrash(this)
            for it = 1:length(this.trash_)
               delete(this.trash_{it}); 
            end
            this.trash_ = {};
        end
        function fqfp = ensureAllRnumbers(~, fqfp, r)
            fqfp = strrep(fqfp, 'r0', sprintf('r%i', r));
        end
        function fqfp = ensureLastRnumber(~, fqfp, r)
            %% ENSURELASTRNUMBER
            %  @param fqfp
            %  @param r integer
            %  @return             ${fqfp}r${r}                     if not fqfp has r[0-9]
            %  @return ${fqfp_upto_r[0-9]}r${r}${fqfp_after_r[0-9]} if fqfp has r[0-9]
            
            assert(isnumeric(r));
            
            startIdx = regexp(fqfp, 'r\d');
            if (~isempty(startIdx))
                fqfp(startIdx(end)+1) = num2str(r);
                return
            end
            fqfp = sprintf('%sr%i', fqfp, r);
        end
        function fqfp = ensureLastRnumbers(~, fqfp, varargin)
            %% ENSURELASTRNUMBERS
            %  @param fqfp
            %  @param ra integer
            %  @param [rb,rc, ...] integer
            %  @return             ${fqfp}r${ra}[r${rb}r${rc} ...]                     if not fqfp has r[0-9]
            %  @return ${fqfp_upto_r[0-9]}r${ra}[r${rb}r${rc} ...]${fqfp_after_r[0-9]} if fqfp has r[0-9]
            
            ip = inputParser;
            addRequired(ip, 'rcell', @(x) iscell(x) && ~isempty(x) && isnumeric(x{1}));
            parse(ip, varargin);
            rcell = ip.Results.rcell;

            concatRs = '';
            for rIdx = 1:length(rcell)
                concatRs = [concatRs sprintf('r%i', rcell{rIdx})]; %#ok<AGROW>
            end
            startIdx = regexp(fqfp, 'r\d');
            if (~isempty(startIdx))                
                fqfp = [fqfp(1:startIdx(end)-1) concatRs fqfp(startIdx(end)+2:end)];
                return
            end
            fqfp = [fqfp concatRs ];
        end
        function fp   = markLastRevisionMarking(~, fp)
            pos = regexp(fp, 'r\d$');
            if (~isempty(pos))
                pos = pos(end);
                fp(pos:pos+1) = 'r0';
            end
        end
        function this = pushTrash(this, t)
            this.trash_ = [this.trash_ t];
        end
        function str  = resolveDest(~, str)
            if (length(str) < 4)
                return
            end
            if (strcmp('op_', str(1:3)))
                str = str(4:end);
            end
        end
        function fqfp = sumTimesActmapf4dfp(this, varargin)
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'fqfpDyn');
            addParameter(ip, 'tag', 'sumt', @ischar);
            parse(ip, varargin{:});
            
            fqfp = sprintf('%s_%s', ip.Results.fqfpDyn, ip.Results.tag);
            if (~lexist_4dfp(fqfp))
                this.buildVisitor.actmapf_4dfp( ...
                    sprintf('%i+', this.imageComposite.length), ip.Results.fqfpDyn, 'options', ['-a' ip.Results.tag]);
            end
        end
        function fqfp = sumTimesImagingContext(this, varargin)
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'fqfpDyn');
            addParameter(ip, 'indicesLogical', this.indicesLogical);          
            parse(ip, varargin{:});
            fqfp = sprintf('%s_sumt%i', ip.Results.fqfpDyn, sum(ip.Results.indicesLogical));
            if (lexist_4dfp(fqfp))
                return
            end
            
            nii = mlfourd.NIfTId([ip.Results.fqfpDyn '.4dfp.ifh']);
            size = nii.size;
            assert(size(4) == length(this.indicesLogical));
            img = zeros(size(1:3));
            for f = 1:length(ip.Results.indicesLogical)
                if (ip.Results.indicesLogical(f))
                    img = img + squeeze(nii.img(:,:,:,f));
                end
            end
            nii.img = img;
            
            nii.saveas([fqfp '.4dfp.ifh']);
        end
        function this = t4imgAll(this, ipr, tag)
            if (this.skipT4imgAll)
                return
            end
            tag = mybasename(tag);
            for f = 1:length(this.indicesLogical)
                if (this.indicesLogical(f))
                    destFp = ipr.dest{f};
                    try
                        this.buildVisitor.t4img_4dfp( ...
                            sprintf('%s_to_%s_t4', destFp, tag), ...
                            destFp, ...
                            'out', sprintf('%s_%s', destFp, tag), ...
                            'options', ['-O' ipr.dest{this.indexOfReference}]);
                    catch ME
                        handwarning(ME);
                    end
                end
            end
        end  
        function fqfp = zeroSlicesOnBoundaries(this, fqfp, n)
            N     = this.imageComposite.length;
            fqfp1 = sprintf('%s_z%ito%i', fqfp,  1, n);
            fqfp2 = sprintf('%s_z%ito%i', fqfp1, N-n+1, N);
            this.buildVisitor.zero_slice_4dfp(fqfp,  'z', 1, n,     fqfp1);
            this.buildVisitor.zero_slice_4dfp(fqfp1, 'z', N-n+1, N, fqfp2);
            fqfp  = fqfp2;            
        end
    end
    
    %% HIDDEN
    
    properties (Hidden)        
        ipResults_ % caching
    end
    
    %% HIDDEN @deprecated
    
    methods (Hidden) 
        function out = t4img_4dfp_0(this, varargin)
            switch (this.NRevisions)
                case 1
                    out = this.t4img_4dfpr1_0(varargin{:});
                case 2
                    out = this.t4img_4dfpr2_0(varargin{:});
                otherwise
                    error('mlfourdfp:unsupportedSwitchCase', ...
                          'CompositeT4ResolveBuilder.t4img_4dfp.NRevisions->%i', this.NRevisions);
            end
        end
        function out = t4img_4dfpr1_0(this, varargin)
            ip = inputParser;
            addRequired(ip, 'src0', @(x) lstrfind(x, this.theImages));
            addRequired(ip, 'src',  @lexist_4dfp);
            addParameter(ip, 'ref', this.theImages{1}, @lexist_4dfp);
            addParameter(ip, 'out', this.buildVisitor.fileprefixT4img(varargin{1:2}), ...
                                            @ischar);
            addParameter(ip, 'options', '', @ischar);
            parse(ip, varargin{:});            
            src0 = ip.Results.src0;
            src  = ip.Results.src;
            ref  = ip.Results.ref;
            opts = ip.Results.options;
            ref_ = [ref 'r1_' this.resolveTag];
            if (~this.buildVisitor.lexist_4dfp(ref_)) %% KLUDGE
                ref_ = ref;
            end
            
            t4  = sprintf('%sr1_to_%s_t4', src0, this.resolveTag);
            in  = src;
            out = sprintf('%sr1_%s', this.clipLastRevisionMarking(src), this.resolveTag);            
            this.buildVisitor.t4img_4dfp(t4, in, 'out', out, 'options', ['-O' ref_ ' ' opts]);
        end
        function out = t4img_4dfpr2_0(this, varargin)
            ip = inputParser;
            addRequired(ip, 'src0', @(x) lstrfind(x, this.theImages));
            addRequired(ip, 'src',  @lexist_4dfp);
            addParameter(ip, 'ref', this.theImages{1}, @lexist_4dfp);
            addParameter(ip, 'out', this.buildVisitor.fileprefixT4img(varargin{1:2}), ...
                                            @ischar);
            addParameter(ip, 'options', '', @ischar);
            parse(ip, varargin{:});            
            src0 = ip.Results.src0;
            src  = ip.Results.src;
            ref  = ip.Results.ref;
            opts = ip.Results.options;
            ref_ = [ref 'r1_' this.resolveTag];
            if (~this.buildVisitor.lexist_4dfp(ref_)) %% KLUDGE
                ref_ = ref;
            end
            
            t4  = sprintf('%sr1_to_%s_t4', src0, this.resolveTag);
            in  = src;
            out = sprintf('%sr1_%s', this.clipLastRevisionMarking(src), this.resolveTag);
            this.buildVisitor.t4img_4dfp(t4, in, 'out', out, 'options', ['-O' ref_ ' ' opts]);
            
            t4  = sprintf('%sr2_to_%s_t4', src0, this.resolveTag);
            in  = sprintf('%sr1_%s', this.clipLastRevisionMarking(src), this.resolveTag);
            out = sprintf('%sr2_%s', this.clipLastRevisionMarking(src), this.resolveTag);            
            this.buildVisitor.t4img_4dfp(t4, in, 'out', out, 'options', ['-O' ref_ ' ' opts]);
        end             
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

