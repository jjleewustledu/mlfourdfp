classdef UmapResolveBuilder < mlfourdfp.AbstractUmapResolveBuilder0
	%% UMAPRESOLVEBUILDER  

	%  $Revision$
 	%  was created 19-Jul-2016 20:08:10
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.
 	
    
    methods (Static)
        function jobs = serialize(c)
            
            assert(isa(c, 'parallel.cluster.Generic'));
            
            hyglys = { 'HYGLY05' 'HYGLY08' 'HYGLY11' 'HYGLY24' 'HYGLY25' };
            jobs   = cell(1, length(hyglys));
            jjlee  = '/scratch/jjlee/raichle/PPGdata/jjlee';
            eSessFqdns = cellfun(@(x) fullfile(jjlee, x), hyglys, 'UniformOutput', false);
            
            for iSess = 1:length(eSessFqdns)
                try
                    jobs{iSess} = c.batch(@mlfourdfp.UmapResolveBuilder.serialBuildUmaps2, 0, {eSessFqdns{iSess}});
                catch ME
                    handwarning(ME);
                end
            end
        end
        function parWilliam
            import mlsystem.*;
            cd(mlraichle.RaichleRegistry.instance.subjectsDir);
            dt = DirTool('*');
            dtFqdns = dt.fqdns;
            parfor idt = 1:length(dtFqdns)
                for v = 1:2
                    cd(dtFqdns{idt});
                    mlfourdfp.UmapResolveBuilder.serialBuildUmaps(pwd, 'iVisit', v);
                end
            end
        end
        function serialBuildUmaps(varargin)
            
            import mlfourdfp.*;
            setenv('PRINTV', '');

            ip = inputParser;
            addRequired( ip, 'sessPath', @isdir);
            addParameter(ip, 'iVisit', 1, @isnumeric);
            addParameter(ip, 'tag', '', @ischar);
            parse(ip, varargin{:});
            cd(ip.Results.sessPath);
            UmapResolveBuilder.diaryv('serialBuildUmaps');
            UmapResolveBuilder.printv('serialBuildUmaps.ip.Results.sessPath->%s\n', ip.Results.sessPath);

            eVisit = mlsystem.DirTool(ip.Results.sessPath);
            if (mlfourdfp.T4ResolveUtilities.isVisit(eVisit.fqdns{ip.Results.iVisit}))
                eTracer = mlsystem.DirTool(eVisit.fqdns{ip.Results.iVisit});
                for iTracer = 1:length(eTracer.fqdns)
                    pth = eTracer.fqdns{iTracer};
                    UmapResolveBuilder.printv('serialBuildUmaps.pth:  %s\n', pth);
                    if ( mlfourdfp.T4ResolveUtilities.isTracer(pth) && ...
                         mlfourdfp.T4ResolveUtilities.isNAC(pth) && ...
                        ~mlfourdfp.T4ResolveUtilities.isEmpty(pth) && ...
                         mlfourdfp.T4ResolveUtilities.matchesTag(ip.Results.sessPath, ip.Results.tag))

                        try
                            cd(pth);
                            UmapResolveBuilder.printv('serialBuildUmaps:  inner try pwd->%s\n', pwd);
                            sessd = mlraichle.SessionData( ...
                                'studyData',   mlraichle.StudyData, ...
                                'sessionPath', ip.Results.sessPath, ...
                                'snumber',     mlfourdfp.T4ResolveUtilities.scanNumber(eTracer.dns{iTracer}), ...
                                'tracer',      mlfourdfp.T4ResolveUtilities.tracerPrefix(eTracer.dns{iTracer}), ...
                                'vnumber',     mlfourdfp.T4ResolveUtilities.visitNumber(eVisit.dns{ip.Results.iVisit}));
                            this = UmapResolveBuilder('sessionData', sessd);
                            this = this.buildUmaps; %#ok<NASGU>                                    
                            save(sprintf('mlfourdfp_UmapResolveBuilder_serialBuildUmaps_this_%s.mat', datestr(now, 30)), 'this');
                        catch ME
                            handwarning(ME);
                        end
                    end
                end
            end
        end
        function reconvertUmapsToE7Format(varargin)
            ip = inputParser;
            addOptional(ip, 'tag', '', @ischar);
            parse(ip, varargin{:});
            
            mlfourdfp.UmapResolveBuilder.triggering('reconvertUmapsToE7Format', 'tag', ip.Results.tag);
        end        
        function viewUmapsOnNative(varargin)
            ip = inputParser;
            addParameter(ip, 'tag', '', @ischar);
            parse(ip, varargin{:});

            import mlfourdfp.* mlsystem.*;
            setenv('PRINTV', '1');
            studyd = mlraichle.StudyData;
            fv = FourdfpVisitor;
            
            eSess = DirTool(fullfile(studyd.subjectsDir, 'HYGLY*'));
            for iSess = 1:length(eSess.fqdns)

                eVisit = DirTool(fullfile(eSess.fqdns{iSess}, 'V*'));
                for iVisit = 1:length(eVisit.fqdns)

                    try
                        pth = fullfile(eVisit.fqdns{iVisit}, sprintf('FDG_%s-NAC', eVisit.dns{iVisit}));
                        pwd0 = pushd(pth);
                        UmapResolveBuilder.printv('viewUmapsOnDynamic:  try pwd->%s\n', pwd);
                        sessd = mlraichle.SessionData( ...
                            'studyData',   studyd, ...
                            'sessionPath', eSess.fqdns{iSess}, ...
                            'ac',          false, ...
                            'tracer',      'FDG', ...
                            'rnumber',     1, ...
                            'vnumber',     T4ResolveUtilities.visitNumber(eVisit.dns{iVisit}));    
                        if (~lexist([sessd.tracerRevision('typ', 'fp') '_b55.4dfp.img']))
                            fv.imgblur_4dfp(sessd.tracerRevision('typ', 'fp'), 5.5);
                        end
                        mlbash(sprintf('fslview %s_b55.4dfp.img -l Cool -b 0,2000 umapsForNAC.4dfp.img -t 0.7 -b 0,0.2', sessd.tracerRevision('typ', 'fp')));
                        mlbash(sprintf('fslview umapSynthv1.4dfp.img %s_sumt.4dfp.img -t 0.5 -l render3', sessd.tracerResolved('typ', 'fp')));
                        popd(pwd0);
                    catch ME
                        handwarning(ME);
                    end
                end
            end
        end 
    end
    
	methods 
		  
 		function this = UmapResolveBuilder(varargin)
 			%% UMAPRESOLVEBUILDER
 			%  Usage:  this = UmapResolveBuilder()

 			this = this@mlfourdfp.AbstractUmapResolveBuilder0(varargin{:});
            this.finished = mlpipeline.Finished(this, 'path', this.logPath, 'tag', lower(this.sessionData.tracer));
        end        
        
        function  this     = buildUmaps(this, varargin)
            this = this.buildF18Umaps(varargin{:});
        end
        function  this     = buildF18Umaps(this, varargin)
            ip = inputParser;
            addParameter(ip, 'indicesLogical', this.indicesLogical_, @isnumeric);
            parse(ip, varargin{:});
            this.indicesLogical_ = ip.Results.indicesLogical;
            
            cd(this.sessionData.fdgNACLocation);
            if (this.isfinished)
                return
            end
            this.convertUmapTo4dfp;
            %this.standardizeResolvedNames;
            
                  ctm    = this.buildCTMasked;
            [this,ctm]   = this.alignCTToSumtResolved(ctm); 
                  ctm    = this.rescaleCT(ctm);
                  umap   = this.buildCarneyUmap(ctm); %%%umap = this.sessionData.umapSynth('typ','fqfp'); % for incremental testing
            [this,umaps] = this.alignUmapToNACFrames(umap);            
            this.pasteImages(umaps, this.umapsForNAC);
            this.convertUmapsToE7Format(umaps);
            this.teardownBuildUmaps;
        end
        function  fqfp     = standardizeResolvedNames(this)
            if (~lexist(this.sessionData.fdgNACResolved('typ', 'fqfn', 'rnumber', this.NRevisions)) && ...
                 lexist(this.sessionData.fdgNACResolved0('fqfn', 'indexMin', this.indexMin, 'indexMax', this.indexMax, 'rnumber', this.NRevisions)))
                this.buildVisitor.copy_4dfp( ...
                    this.sessionData.fdgNACResolved0('fqfp', 'indexMin', this.indexMin, 'indexMax', this.indexMax, 'rnumber', this.NRevisions), ...
                    this.sessionData.fdgNACResolved('typ', 'fqfp', 'rnumber', this.NRevisions));
            end
            fqfp = this.sessionData.fdgNACResolved('typ', 'fqfp', 'rnumber', this.NRevisions);
        end
        function [this,ct] = alignCTToSumtResolved(this, ct)
            sumtResolved = this.sumTimes(this.sessionData.fdgNACResolved('typ', 'fqfp'));
            mpr          = this.mprOnPetSumt(this.sessionData.mpr('typ', 'fqfp'), sumtResolved);   
            ct           = this.ctOnPetSumt(ct, sumtResolved);            
            ct           = this.resolveSequence(sumtResolved, mpr, ct);
        end
        function  umap     = buildCarneyUmap(this, varargin)
            %% BUILDSYNTHETICUMAP follows Carney, et al. Med. Phys. 33(4) 2006 976-983.
            %  @param ct   is the (fully-qualified) fileprefix of the rescaled CT.
            %  @param umap is the (fully-qualified) fileprefix of the product umap.
            %  @returns umap product as f.q. fileprefix.
            
            import mlfourdfp.*;
            ip = inputParser;
            addOptional(ip, 'rescaledCT', this.buildCTMasked, @FourdfpVisitor.lexist_4dfp);
            addOptional(ip, 'umap', this.sessionData.umapSynth('typ', 'fqfp'), @ischar);
            parse(ip, varargin{:});
            
            umap = ip.Results.umap;
            if (FourdfpVisitor.lexist_4dfp(umap))
                if (this.reuseCarneyUmap)
                    return
                else
                    delete([umap '.4dfp.*']);
                end
            end
            ctic = this.CarneyImagingContext(ip.Results.rescaledCT);
            ctic.saveas([umap '.4dfp.hdr']);
        end
        function             teardownBuildUmaps(this)
            this.teardownLogs;
            this.teardownT4s;
            
            %ensuredir(this.onAtlasPath);
            %movefiles(sprintf('*%s*', this.atlas('typ', 'fp')), this.onAtlasPath);
            %ensuredir(this.resolveSequenceLocation);
            %movefiles(sprintf('%s*', this.resolveSequenceTag), this.resolveSequenceLocation);
            
            %delete([this.resolveSequenceTag '*_frame*4dfp*']);
            %delete(sprintf('%s_on_*.4dfp.*', this.sessionData.ct('typ', 'fp')));  
            %delete(sprintf('%s_on_*.4dfp.*', this.sessionData.ctMasked('typ', 'fp')));
            this.finished.touchFinishedMarker;
        end
        
        %% UTILITY
          
        function fp    = resolveSequenceTag(this)
            fp = sprintf('resolveSequencev%i', this.sessionData.vnumber);
        end
        function umaps = umapsForNAC(~)
            umaps = 'umapsForNAC';
        end
    end
    
    %% PROTECTED
    
    methods (Access = protected)
        function umap = CarneyImagingContext(this, varargin)
            %% CARNEYIMAGINGCONTEXT follows Carney, et al. Med. Phys. 33(4) 2006 976-983.
            %  @param ct   is the (fully-qualified) fileprefix of the rescaled CT.
            %  @returns umap ImagingContext.
            
            import mlfourdfp.*;
            ip = inputParser;
            addOptional(ip, 'ctRescaled', this.buildCTMasked, @FourdfpVisitor.lexist_4dfp);
            %addOptional(ip, 'ctMask', this.sessionData.ctMask('typ', '4dfp.ifh'), @FourdfpVisitor.lexist_4dfp);
            parse(ip, varargin{:});
            
            ct  = mlfourd.ImagingContext([ip.Results.ctRescaled '.4dfp.hdr']);
            ct  = ct.numericalNiftid;
            %ctm = mlfourd.ImagingContext([ip.Results.ctMask '.4dfp.hdr']);
            %ctm = ctm.numericalNiftid;
            
            lowHU    = ct.uthresh(this.CarneyBP); 
            lowMask  = lowHU.binarized;
            highHU   = ct.thresh(this.CarneyBP);
            highMask = highHU.binarized;
            
            lowHU = (lowHU + 1000) * 9.6e-5;
            lowHU =  lowHU .* lowMask;
            
            highHU = (highHU + 1000) * this.CarneyA + this.CarneyB;
            highHU =  highHU .* highMask;

            umap = lowHU + highHU;
            umap = umap .* (umap > 0);
        end   
        function a    = CarneyA(this)
            switch (this.ct_kVp)
                case 80
                    a = 3.64e-5;
                case 100
                    a = 4.43e-5;
                case 110
                    a = 4.92e-5;
                case 120
                    a = 5.10e-5;
                case 130
                    a = 5.51e-5;
                case 140
                    a = 5.64e-5;
                otherwise
                    error('mlfourdfp:valueOutOfBounds', 'UMapResolveBuilder.ct_kVp->%g is not supported', this.ct_kVp);
            end
        end
        function b    = CarneyB(this)
            switch (this.ct_kVp)
                case 80
                    b = 6.26e-2;
                case 100
                    b = 5.44e-2;
                case 110
                    b = 4.88e-2;
                case 120
                    b = 4.71e-2;
                case 130
                    b = 4.24e-2;
                case 140
                    b = 4.08e-2;
                otherwise
                    error('mlfourdfp:valueOutOfBounds', 'UMapResolveBuilder.ct_kVp->%g is not supported', this.ct_kVp);
            end
        end
        function bp   = CarneyBP(this)
            switch (this.ct_kVp)
                case 80
                    bp = 50;
                case 100
                    bp = 52;
                case 110
                    bp = 43;
                case 120
                    bp = 47;
                case 130
                    bp = 37;
                case 140
                    bp = 30;
                otherwise
                    error('mlfourdfp:valueOutOfBounds', 'UMapResolveBuilder.ct_kVp->%g is not supported', this.ct_kVp);
            end
        end
        function ct   = resolveSequence(this, sumtResolved, mpr, ct)
            pasted = fullfile(myfileparts(ct), this.resolveSequenceTag);
            images = {sumtResolved ct mpr};
            assert(lstrfind(images{this.mprFromFdgNacResolveSequence}, this.sessionData.mpr('typ', 'fp')));
            this.pasteImages(images, pasted);
            this = this.resolve( ...
                'dest', [pasted '_resolved'], ...
                'source', pasted, ...
                'sourceBlur', [this.blurArg this.blurArg this.blurArg], ...
                'indicesLogical', [1 1 1], ...
                'NRevisions', this.resolveSequenceNRevisions);
            resolved = sprintf('%sr%i_resolved', pasted, this.NRevisions);
            this.buildVisitor.extract_frame_4dfp(resolved, 2);
            this.buildVisitor.delete_4dfp(ct);
            this.buildVisitor.move_4dfp([resolved '_frame2'], ct);
        end   
    end
    
    %% PRIVATE
    
    methods (Static, Access = private)
        function this = triggering(varargin)
            
            studyd = mlraichle.StudyData;
            
            ip = inputParser;
            addRequired( ip, 'methodName', @ischar);
            addParameter(ip, 'subjectsDir', studyd.subjectsDir, @isdir);
            addParameter(ip, 'tag', '', @ischar);
            parse(ip, varargin{:});
            if (~strcmp(ip.Results.subjectsDir, studyd.subjectsDir))
                studyd.subjectsDir = ip.Results.subjectsDir;
            end
            
            import mlsystem.* mlfourdfp.* ;
            UmapResolveBuilder.printv('triggering.ip.Results:  %s\n', struct2str(ip.Results));            
            eSess = DirTool(ip.Results.subjectsDir);
            UmapResolveBuilder.printv('triggering.eSess:  %s\n', cell2str(eSess.dns));
            for iSess = 1:length(eSess.fqdns)

                eVisit = DirTool(eSess.fqdns{iSess});
                UmapResolveBuilder.printv('triggering.eVisit:  %s\n', cell2str(eVisit.dns));
                for iVisit = 1:length(eVisit.fqdns)
                    
                        eTracer = DirTool(eVisit.fqdns{iVisit});
                        UmapResolveBuilder.printv('triggering.eTracer:  %s\n', cell2str(eTracer.dns));
                        for iTracer = 1:length(eTracer.fqdns)

                            pth___ = eTracer.fqdns{iTracer};
                            UmapResolveBuilder.printv('triggering.pth:  %s\n', pth___);
                            if (T4ResolveUtilities.matchesTag(pth___, ip.Results.tag) && ...
                                T4ResolveUtilities.isVisit(pth___) && ...
                                T4ResolveUtilities.isTracer(pth___) && ...
                                T4ResolveUtilities.isNAC(pth___) && ...
                               ~T4ResolveUtilities.isEmpty(pth___))
                                try
                                    sessd = mlraichle.SessionData( ...
                                        'studyData',   studyd, ...
                                        'sessionPath', eSess.fqdns{iSess}, ...
                                        'snumber',     T4ResolveUtilities.scanNumber(eTracer.dns{iTracer}), ...
                                        'tracer',      T4ResolveUtilities.tracerPrefix(eTracer.dns{iTracer}), ...
                                        'vnumber',     T4ResolveUtilities.visitNumber(eVisit.dns{iVisit}));
                                    this = UmapResolveBuilder('sessionData', sessd);
                                    this.(ip.Results.methodName);   
                                catch ME
                                    handwarning(ME);
                                end
                            end
                        end
                end                
            end            
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

