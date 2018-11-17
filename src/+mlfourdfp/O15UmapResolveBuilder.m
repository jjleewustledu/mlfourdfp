classdef O15UmapResolveBuilder < mlfourdfp.AbstractUmapResolveBuilder
	%% O15UMAPRESOLVEBUILDER  

	%  $Revision$
 	%  was created 24-Oct-2016 21:03:18
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.
 	

    methods (Static)
        function prepareWilliam
            import mlsystem.*;
            cd(mlraichle.RaichleRegistry.instance.subjectsDir);            
            dt = DirTool('*');
            dtFqdns = dt.fqdns;            
            studyd = mlraichle.StudyData;
            parfor v = 1:2
                for idt = 1:length(dtFqdns)
                    cd(dtFqdns{idt});
                    sessd = mlraichle.SessionData('studyData', studyd, 'sessionPath', pwd);
                    sessd.vnumber = v;
                    o15urb = mlfourdfp.O15UmapResolveBuilder('sessionData', sessd);
                    o15urb = o15urb.loadSessionDataCache({'HO' 'OO' 'OC'});
                    for sdc = 1:length(o15urb.sessionDataCache)
                        try
                            o15urb.buildTracerNAC(o15urb.sessionDataCache{sdc});
                        catch ME
                            handwarning(ME);
                        end
                    end
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
                    mlfourdfp.O15UmapResolveBuilder.buildAllO15Umaps(pwd, 'iVisit', v);
                end
            end
        end
        function reconvertUmapToE7Format(varargin)
            ip = inputParser;
            addOptional(ip, 'tag', '', @ischar);
            parse(ip, varargin{:});
            
            mlfourdfp.O15UmapResolveBuilder.triggering('convertUmapToE7Format', 'tag', ip.Results.tag);
        end   
        function buildAllO15Umaps(varargin)
            
            import mlfourdfp.*;
            setenv('PRINTV', '');

            ip = inputParser;
            addRequired( ip, 'sessPath', @isdir);
            addParameter(ip, 'iVisit', 1, @isnumeric);
            addParameter(ip, 'tag', '', @ischar);
            parse(ip, varargin{:});
            cd(ip.Results.sessPath);
            O15UmapResolveBuilder.diaryv('serialBuildUmaps');
            O15UmapResolveBuilder.printv('serialBuildUmaps.ip.Results.sessPath->%s\n', ip.Results.sessPath);

            eVisit = mlsystem.DirTool(ip.Results.sessPath);
            if (mlfourdfp.T4ResolveUtilities.isVisit(eVisit.fqdns{ip.Results.iVisit}))
                
                pth = eVisit.fqdns{ip.Results.iVisit};
                O15UmapResolveBuilder.printv('serialBuildUmaps.pth:  %s\n', pth);

                if (~mlfourdfp.T4ResolveUtilities.isEmpty(pth) && ...
                     mlfourdfp.T4ResolveUtilities.matchesTag(ip.Results.sessPath, ip.Results.tag))

                    try
                        cd(pth);
                        O15UmapResolveBuilder.printv('serialBuildUmaps:  inner try pwd->%s\n', pwd);
                        sessd = mlraichle.SessionData( ...
                            'studyData',   mlraichle.StudyData, ...
                            'sessionPath', ip.Results.sessPath, ...
                            'vnumber',     mlfourdfp.T4ResolveUtilities.visitNumber(eVisit.dns{ip.Results.iVisit}));
                        this = O15UmapResolveBuilder('sessionData', sessd);
                        this = this.buildUmaps; %#ok<NASGU>                                    
                        save(sprintf('mlfourdfp_O15UmapResolveBuilder_buildAllO15Umaps_this_%s.mat', datestr(now, 30)), 'this');
                    catch ME
                        handwarning(ME);
                    end
                end
            end
        end 
        function buildAllO15AfterAC
            import mlfourdfp.*;
            O15UmapResolveBuilder.triggeringSessionDataCache('buildO15AfterAC', 'tracerExpr', {'HO' 'OO'});
            O15UmapResolveBuilder.triggeringSessionDataCache('buildO15AfterAC', 'tracerExpr', {'OC'});            
        end
        function prepareViewAllAfterAC
            mlfourdfp.O15UmapResolveBuilder.triggeringSessionDataCache('prepareViewAfterAC', 'tracerExpr', {'HO' 'OO' 'OC'});
        end
        function t4imgAllAfterAC
            mlfourdfp.O15UmapResolveBuilder.triggeringSessionDataCache('t4imgAfterAC', 'tracerExpr', {'HO' 'OO' 'OC'});
        end
        function viewAllAfterAC
            mlfourdfp.O15UmapResolveBuilder.triggeringVisits('viewAfterAC');
        end
        function viewAllUmapOnNAC(varargin)
            ip = inputParser;
            addOptional(ip, 'tag', '', @ischar);
            parse(ip, varargin{:});
            
            mlfourdfp.O15UmapResolveBuilder.triggering('viewUmapOnNAC', 'tag', ip.Results.tag);
        end   
        function viewAllUmapMprEarlyAC(varargin)
            ip = inputParser;
            addOptional(ip, 'tag', '', @ischar);
            parse(ip, varargin{:});
            
            mlfourdfp.O15UmapResolveBuilder.triggering('viewUmapMprEarlyAC', 'tag', ip.Results.tag);
        end
    end
    
	methods
 		function this = O15UmapResolveBuilder(varargin)
 			%% O15UMAPRESOLVEBUILDER
 			%  Usage:  this = O15UmapResolveBuilder()

 			this = this@mlfourdfp.AbstractUmapResolveBuilder(varargin{:});
            this.f18UmapResolveBuilder_ = mlfourdfp.CarneyUmapBuilder(varargin{:});
            sessd = this.sessionData;
            sessd.tracer = 'FDG';
            assert(this.f18UmapResolveBuilder_.isfinished);
            cd(sessd.vLocation);
            this.finished_ = mlpipeline.Finished(this, 'path', this.logPath, 'tag', lower(this.sessionData.tracer));
        end
        
        function tf = isfinished(~, sessd)
            if (mlfourdfp.O15UmapResolveBuilder.REPLACE_COMPLETED)
                tf = false;
                return
            end
            assert(isa(sessd, 'mlpipeline.SessionData'));
            
            sessd_oo = sessd;
            sessd_oo.tracer = 'OO';
            sessd_oo.snumber = 2;
            ooUrb = mlfourdfp.O15UmapResolveBuilder('sessionData', sessd_oo);
            
            sessd_oc = sessd;
            sessd_oc.tracer = 'OC';
            sessd_oc.snumber = 2;
            ocUrb = mlfourdfp.O15UmapResolveBuilder('sessionData', sessd_oc);
            
            tf = lexist(ooUrb.finished.finishedMarkerFilename, 'file') && lexist(ocUrb.finished.finishedMarkerFilename, 'file');
        end
        function tf = isfinishedOC(~, sessd)
            if (mlfourdfp.O15UmapResolveBuilder.REPLACE_COMPLETED)
                tf = false;
                return
            end
            assert(isa(sessd, 'mlpipeline.SessionData'));
            sessd.tracer = 'OC';
            sessd.snumber = 2;
            this = mlfourdfp.O15UmapResolveBuilder('sessionData', sessd);
            tf = lexist(this.finished.finishedMarkerFilename, 'file');
        end
        function tf = isfinishedOO(~, sessd)
            if (mlfourdfp.O15UmapResolveBuilder.REPLACE_COMPLETED)
                tf = false;
                return
            end
            assert(isa(sessd, 'mlpipeline.SessionData'));
            sessd.tracer = 'OO';
            sessd.snumber = 2;
            this = mlfourdfp.O15UmapResolveBuilder('sessionData', sessd);
            tf = lexist(this.finished.finishedMarkerFilename, 'file');
        end
        function this            = buildUmaps(this, varargin)
%             if (~this.isfinishedOC(this.sessionData))
%                 this.blurArg_ = 1.5;
%                 this = this.loadSessionDataCache({'OC'});
%                 this = this.buildO15Umaps(varargin{:});
%             end
            if (~this.isfinishedOO(this.sessionData))
                this = this.loadSessionDataCache({'HO' 'OO'});
                this = this.buildO15Umaps(varargin{:});
            end
        end     
        function [this,resolved] = buildO15Umaps(this, varargin)
            o15Nacs = {};
            for sdc = 1:length(this.sessionDataCache)
                this.sessionData_ = this.sessionDataCache{sdc};
                o15Nacs = [o15Nacs {this.buildTracerNAC}]; %#ok<AGROW>
            end
            resolved =  this.alignO15NACs(o15Nacs);             
            [this,umaps] = this.alignUmapToNACs(o15Nacs);
            this.reconstituteImages(struct('fqfps', umaps,   'dest1', this.umapsForNAC));
            this.reconstituteImages(struct('fqfps', o15Nacs, 'dest1', this.tracerForNAC));
            this.convertUmapsToE7Format(umaps);
            this.teardownBuildUmaps;
        end
        function                   teardownBuildUmaps(this)
            this.teardownLogs;
            this.teardownT4s;
            
            %delete([this.resolveSequenceTag '*_frame*4dfp*']);
            %delete([this.resolveSequenceTag '*_frame*nii.gz']);
            %ensuredir(this.resolveSequenceLocation);
            %movefiles(sprintf('%s*', this.resolveSequenceTag), this.resolveSequenceLocation);
            this.finished.touchFinishedMarker;
        end 
        function this            = buildO15AfterAC(this, sdc)
            assert(isnumeric(sdc));
            
            pthAC = this.sessionData.tracerLocation;
            if (isdir(pthAC))
                rmdir(pthAC, 's');
            end
            ensuredir(pthAC);
            fprintf('O15UmapResolveBuilder.buildO15AfterAC:  working with -> %s\n', pthAC);
                               
            trLM = this.sessionData.tracerListmodeSif('typ', 'fp');
            trRev = this.sessionData.tracerRevision('typ', 'fp');
            bv = this.buildVisitor;
            %if (~lexist(fullfile(pthAC, [trRev '.4dfp.hdr']), 'file'))
                pwd0 = pushd(this.sessionData.tracerListmodeLocation); 
                bv.sif_4dfp(trLM);
                bv.cropfrac_4dfp(0.5, trLM, trRev);
                bv.move_4dfp(trRev, fullfile(pthAC, trRev));
                delete('*.4dfp.*')
                popd(pwd0);
            %end

            pwd1 = pushd(pthAC);
            Ncache = length(this.sessionDataCache);
            bv.lns(this.tracerToResolvedT4(sdc+1));
            bv.lns(this.resolvedToTracerT4(sdc+1));
            bv.lns(this.tracerToMprT4(sdc+1, Ncache+3));
            bv.lns(this.mprToTracerT4(Ncache+3, sdc+1));
%            bv.imgblur_all_4dfp(5.5);
            popd(pwd1);
        end
        function this            = t4imgAfterAC(this, ~)
            pthAC = this.sessionData.tracerLocation;
            trRev = this.sessionData.tracerRevision('typ', 'fp');
            pwd0 = pushd(pthAC);
            fprintf('O15UmapResolveBuilder.t4imgAfterAC:  working in -> %s\n', fullfile(pthAC, trRev));
                               
            t4 = sprintf('%s_to_resolved_t4', trRev);
            resolved = this.sessionData.fdgNACResolved('typ', 'fqfp');
            this.buildVisitor.t4img_4dfp(t4, trRev, [trRev '_on_resolved'], 'options', ['-O' resolved]);
            popd(pwd0);
        end
        function this            = prepareViewAfterAC(this, ~)
            sessd = this.sessionData;
            pthAC = sessd.tracerLocation;
            trRev = sessd.tracerRevision('typ', 'fp');
            pwd0 = pushd(pthAC);
            fprintf('O15UmapResolveBuilder.prepareViewAfterAC:  working with -> %s\n', fullfile(pthAC, trRev));
                    
            t4 = fullfile(pthAC, sprintf('%s_to_resolved_t4', trRev));
            trRevSumt = this.sumTimes(trRev);
            resolvedSumt = [sessd.fdgNACResolved('typ', 'fqfp') '_sumt'];
            this.buildVisitor.t4img_4dfp(t4, trRevSumt, 'out', [trRevSumt '_on_resolved'], 'options', ['-O' resolvedSumt]);
            
            popd(pwd0);
        end
        function this            = viewAfterAC(this)
            sessd = this.sessionData;
            Vn   = sessd.vLocation('typ', 'folder');
            oc1  = fullfile(sessd.vLocation, ['OC1_' Vn '-AC'], sprintf('oc1v%ir1_sumt_on_resolved.4dfp.img', sessd.vnumber));
            oc2  = fullfile(sessd.vLocation, ['OC2_' Vn '-AC'], sprintf('oc2v%ir1_sumt_on_resolved.4dfp.img', sessd.vnumber));
            oo1  = fullfile(sessd.vLocation, ['OO1_' Vn '-AC'], sprintf('oo1v%ir1_sumt_on_resolved.4dfp.img', sessd.vnumber));
            oo2  = fullfile(sessd.vLocation, ['OO2_' Vn '-AC'], sprintf('oo2v%ir1_sumt_on_resolved.4dfp.img', sessd.vnumber));
            ho1  = fullfile(sessd.vLocation, ['HO1_' Vn '-AC'], sprintf('ho1v%ir1_sumt_on_resolved.4dfp.img', sessd.vnumber));
            ho2  = fullfile(sessd.vLocation, ['HO2_' Vn '-AC'], sprintf('ho2v%ir1_sumt_on_resolved.4dfp.img', sessd.vnumber));
            imgs = {};
            frmt = '';
            if (lexist(oc1, 'file')); imgs = [imgs oc1]; frmt = [frmt ' %s -l render3']; end
            if (lexist(oc2, 'file')); imgs = [imgs oc2]; frmt = [frmt ' %s -l render3']; end
            if (lexist(oo1, 'file')); imgs = [imgs oo1]; frmt = [frmt ' %s']; end
            if (lexist(oo2, 'file')); imgs = [imgs oo2]; frmt = [frmt ' %s']; end
            if (lexist(ho1, 'file')); imgs = [imgs ho1]; frmt = [frmt ' %s']; end
            if (lexist(ho2, 'file')); imgs = [imgs ho2]; frmt = [frmt ' %s']; end
            mlbash(sprintf(['fslview ' frmt ' %s.4dfp.img -l Cool -b 0,5000 %s.4dfp.img -t 0.5'], ...
                   imgs{:}, this.fdgOnResolvedSumt, this.mprOnResolved));
        end
        function fdg             = fdgOnResolvedSumt(this)
            pwd0 = pushd(this.sessionData.fdgACRevision('typ', 'path'));
            fdg = [this.sessionData.fdgACRevision('typ', 'fqfp') '_on_resolved'];   
            if (~lexist([fdg '_sumt.4dfp.img']))
                fdg = this.sumTimes(fdg);
            else
                fdg = [fdg '_sumt'];
            end
            popd(pwd0);
        end
        function mpr             = mprOnResolved(this)
            sessd = this.sessionData;
            mpr = [sessd.mpr('typ', 'fqfp') '_on_' sessd.fdgNACResolved('typ', 'fp')];
        end
        function                   viewUmapOnNAC(this)
            sessd = this.sessionData;
            sessd.rnumber = 1;
            umap = fullfile( ...
                sessd.tracerNACLocation, ...
                sprintf('%s_on_%s.4dfp.img', sessd.umapSynth('typ', 'fp'), sessd.tracerNACRevision('typ', 'fp')));
            tracerRevision = sessd.tracerRevision('typ', '.4dfp.img');
            mlbash(sprintf('fslview %s %s -l render3 -t 0.6 -b 0.01,0.13', ...
                tracerRevision, umap));
        end
        function                   viewUmapMprEarlyAC(this)
            sessd = this.sessionData;
            sessd.rnumber = 1;
            sessdAC = sessd;
            sessdNAC = sessd;
            sessdAC.attenuationCorrected = true;
            sessdNAC.attenuationCorrected = false;
            sessdOO = sessdNAC;
            sessdOO.tracer = 'OO';
            sessdOO.snumber = 2;
            if (~isdir(sessdOO.tracerLocation))
                sessdOO.snumber = 1;
                assert(isdir(sessdOO.tracerLocation));
            end
            tr0 = sessdAC.tracerRevision('typ', 'fqfp');
            seq = fullfile( ...
                sessdOO.tracerLocation, ...
                upperFirst(this.resolveSequenceTag), this.resolveSequenceResolved('typ', 'fp'));
            ipr = struct('dest', seq, 'frame', 1);
            fdg = this.imageComposite.lazyExtractImage(ipr);
            ipr.frame = this.imageComposite.readLength(seq);
            mpr = this.imageComposite.lazyExtractImage(ipr);
            tr  = this.framedSummedImaging(tr0, ones(1,20));
            umap = fullfile( ...
                sessd.tracerNACLocation, ...
                sprintf('%s_on_%s', sessdNAC.umapSynth('typ', 'fp'), sessdNAC.tracerNACRevision('typ', 'fp')));
            mlbash(sprintf('fslview %s.4dfp.img %s.4dfp.img %s.4dfp.img %s.4dfp.img -l render3 -t 0.6 -b 0.01,0.13', ...
                fdg, tr, mpr, umap));
        end
        function fqfp = framedImaging(~, fqfp0, indicesLogical)
            assert(isnumeric(indicesLogical) && 1 == max(indicesLogical));
            ic0     = mlfourd.ImagingContext([fqfp0 '.4dfp.hdr']);
            nii0    = ic0.niftid;
            nii     = nii0;
            nii.img = zeros(nii0.size(1), nii0.size(2), nii0.size(3), sum(indicesLogical));
            u       = 0;
            for t = 1:length(indicesLogical)
                u = u + 1;
                if (indicesLogical(t))
                    nii.img(:,:,:,u) = nii0.img(:,:,:,t);
                end
            end
            [~,indexMin] = max(indicesLogical);
            fqfp = sprintf('%s_frames%ito%i', fqfp0, indexMin, sum(indicesLogical));
            nii.noclobber = false;
            nii.saveas([fqfp '.4dfp.hdr']);
        end
        function fqfp = framedSummedImaging(~, fqfp0, indicesLogical)
            assert(isnumeric(indicesLogical) && 1 == max(indicesLogical));
            ic0     = mlfourd.ImagingContext([fqfp0 '.4dfp.hdr']);
            nii0    = ic0.niftid;
            nii     = nii0;
            nii.img = zeros(nii0.size(1), nii0.size(2), nii0.size(3));
            u       = 0;
            for t = 1:length(indicesLogical)
                u = u + 1;
                if (indicesLogical(t))
                    nii.img = nii.img + squeeze(nii0.img(:,:,:,t));
                end
            end
            [~,indexMin] = max(indicesLogical);
            fqfp = sprintf('%s_frames%ito%i_sumt', fqfp0, indexMin, sum(indicesLogical));
            nii.noclobber = false;
            nii.saveas([fqfp '.4dfp.hdr']);
        end
 	end 
    
    methods (Static) %%% DEBUG , Access = private)
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
            O15UmapResolveBuilder.printv('triggering.ip.Results:  %s\n', struct2str(ip.Results));            
            eSess = DirTool(ip.Results.subjectsDir);
            O15UmapResolveBuilder.printv('triggering.eSess:  %s\n', cell2str(eSess.dns));
            for iSess = 1:length(eSess.fqdns)

                eVisit = DirTool(eSess.fqdns{iSess});
                O15UmapResolveBuilder.printv('triggering.eVisit:  %s\n', cell2str(eVisit.dns));
                for iVisit = 1:length(eVisit.fqdns)
                        
                    eTracer = DirTool(eVisit.fqdns{iVisit});
                    O15UmapResolveBuilder.printv('triggering.eTracer:  %s\n', cell2str(eTracer.dns));
                    for iTracer = 1:length(eTracer.fqdns)

                        pth___ = eTracer.fqdns{iTracer};
                        O15UmapResolveBuilder.printv('triggering.pth:  %s\n', pth___);
                        if (T4ResolveUtilities.matchesTag(pth___, ip.Results.tag) && ...
                            T4ResolveUtilities.isVisit(pth___) && ...
                            T4ResolveUtilities.isTracer(pth___, {'HO' 'OO'}) && ...
                            T4ResolveUtilities.isNAC(pth___) && ...
                           ~T4ResolveUtilities.isEmpty(pth___) && ...
                           ~T4ResolveUtilities.isConverted(pth___) && ...
                           ~T4ResolveUtilities.isBackup(pth___))
                            try
                                sessd = mlraichle.SessionData( ...
                                    'studyData',   studyd, ...
                                    'sessionPath', eSess.fqdns{iSess}, ...
                                    'snumber',     T4ResolveUtilities.scanNumber(eTracer.dns{iTracer}), ...
                                    'tracer',      T4ResolveUtilities.tracerPrefix(eTracer.dns{iTracer}), ...
                                    'vnumber',     T4ResolveUtilities.visitNumber(eVisit.dns{iVisit}));
                                this = O15UmapResolveBuilder('sessionData', sessd);
                                this.(ip.Results.methodName);   
                            catch ME
                                handwarning(ME);
                            end
                        end
                    end
                end                
            end            
        end
        function this = triggeringSessionDataCache(varargin)
            
            ip = inputParser;
            addRequired( ip, 'methodName', @ischar);
            addParameter(ip, 'ac', true, @islogical);
            addParameter(ip, 'sessionsExpr', 'HYGLY*', @ischar);
            addParameter(ip, 'studyd', mlraichle.StudyData, @(x) isa(x, 'mlpipeline.IStudyHandle'));
            addParameter(ip, 'tracerExpr', {'HO' 'OO'}, @iscell);
            addParameter(ip, 'visitExpr', {'V1' 'V2'}, @iscell);
            parse(ip, varargin{:});
            
            import mlsystem.* mlfourdfp.*;     
            eSess = DirTool(fullfile(ip.Results.studyd.subjectsDir, ip.Results.sessionsExpr));
            for iSess = 1:length(eSess.fqdns)

                eVisit = DirTool(eSess.fqdns{iSess});
                for iVisit = 1:length(eVisit.fqdns)
                        
                    if (lstrfind(ip.Results.visitExpr, eVisit.dns{iVisit}))

                        pwd0_ = pushd(eVisit.fqdns{iVisit});                        
                        sessd_ = mlraichle.SessionData( ...
                            'studyData',   ip.Results.studyd, ...
                            'sessionPath', eSess.fqdns{iSess}, ...
                            'ac',          ip.Results.ac, ...
                            'vnumber',     T4ResolveUtilities.visitNumber(eVisit.dns{iVisit}));
                        this = O15UmapResolveBuilder('sessionData', sessd_);
                        this = this.loadSessionDataCache(ip.Results.tracerExpr);
                        for sdc = 1:length(this.sessionDataCache)

                            try
                                this.sessionData_ = this.sessionDataCache{sdc};
                                this = this.(ip.Results.methodName)(sdc);
                            catch ME
                                handwarning(ME);
                            end
                        end                        
                        popd(pwd0_)                    
                    end
                end                
            end
        end
        function this = triggeringVisits(varargin)
            
            ip = inputParser;
            addRequired( ip, 'methodName', @ischar);
            addParameter(ip, 'ac', true, @islogical);
            addParameter(ip, 'sessionsExpr', 'HYGLY*', @ischar);
            addParameter(ip, 'studyd', mlraichle.StudyData, @(x) isa(x, 'mlpipeline.IStudyHandle'));
            addParameter(ip, 'visitExpr', {'V1' 'V2'}, @iscell);
            parse(ip, varargin{:});
            
            import mlsystem.* mlfourdfp.*;
            eSess = DirTool(fullfile(ip.Results.studyd.subjectsDir, ip.Results.sessionsExpr));
            for iSess = 1:length(eSess.fqdns)

                eVisit = DirTool(eSess.fqdns{iSess});
                for iVisit = 1:length(eVisit.fqdns)
                        
                    if (lstrfind(ip.Results.visitExpr, eVisit.dns{iVisit}))

                        pwd0_ = pushd(eVisit.fqdns{iVisit});  
                        try
                            sessd_ = mlraichle.SessionData( ...
                                'studyData',   ip.Results.studyd, ...
                                'sessionPath', eSess.fqdns{iSess}, ...
                                'ac',          ip.Results.ac, ...
                                'vnumber',     T4ResolveUtilities.visitNumber(eVisit.dns{iVisit}));
                            this = O15UmapResolveBuilder('sessionData', sessd_);
                            this = this.(ip.Results.methodName);
                        catch ME
                            handwarning(ME);
                        end
                        popd(pwd0_)                        
                    end
                end                
            end
        end
    end
    
    methods
        function t4 = mprToTracerT4(this, frameMpr, frameTr)
            assert(isnumeric(frameTr));
            rsTag = this.f18UmapResolveBuilder_.resolveSequenceTag;
            t4 = this.mapcycle(frameMpr, frameTr, ...
                               sprintf('%sr1_frame%i_to_%s_t4', ...
                                       rsTag, frameMpr, ...
                                       this.sessionData.tracerRevision('typ', 'fp')));
            t4 = this.leftMultiplyMprT4(t4);
        end
        function t4 = tracerToMprT4(this, frameTr, frameMpr)
            assert(isnumeric(frameTr));
            rsTag = this.f18UmapResolveBuilder_.resolveSequenceTag;
            t4 = this.mapcycle(frameTr, frameMpr, ...
                               sprintf('%s_to_%sr1_frame%i_t4', ...
                                       this.sessionData.tracerRevision('typ', 'fp'), ...
                                       rsTag, frameMpr));   
            t4 = this.rightMultiplyMprT4(t4);         
        end
        function t4 = tracerToResolvedT4(this, frameTr)
            assert(isnumeric(frameTr));
            t4 = this.mapcycle(frameTr, 1, ...
                               sprintf('%s_to_resolved_t4', this.sessionData.tracerRevision('typ', 'fp')));            
        end
        function t4 = resolvedToTracerT4(this, frameTr)
            assert(isnumeric(frameTr));
            t4 = this.mapcycle(1, frameTr, ...
                               sprintf('resolved_to_%s_t4', this.sessionData.tracerRevision('typ', 'fp')));            
        end
        function t4 = leftMultiplyMprT4(this, t4)
            mpr_to_fdgNacResolveSequence_t4 = ...
                fullfile(this.sessionData.fdgNACT4Location, ...
                         sprintf('mpr_to_%s_sumt_t4', this.sessionData.fdgNACResolved('typ', 'fp')));
            t4 = this.buildVisitor.t4_mul(mpr_to_fdgNacResolveSequence_t4, t4);
        end
        function t4 = rightMultiplyMprT4(this, t4)
            fdgNacResolveSequence_to_mpr_t4 = ...
                fullfile(this.sessionData.fdgNACT4Location, ...
                         sprintf('%s_sumt_to_mpr_t4', this.sessionData.fdgNACResolved('typ', 'fp')));
            t4 = this.buildVisitor.t4_mul(t4, fdgNacResolveSequence_to_mpr_t4);
        end
    end
    
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

