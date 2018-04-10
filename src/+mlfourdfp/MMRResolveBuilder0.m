classdef MMRResolveBuilder0 < mlfourdfp.T4ResolveBuilder0
	%% MMRRESOLVEBUILDER0  

	%  $Revision$
 	%  was created 01-Nov-2016 19:09:02
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.
 	

	properties (Constant) 		
        CLUSTER_HOSTNAME = 'dtn01.chpc.wustl.edu'
        CLUSTER_SUBJECTS_DIR = '/scratch/jjlee/raichle/PPGdata/jjlee'
    end
    
    properties   
        ct_kVp = 120
        ct_rescaleSlope = 1
        ct_rescaleIntercept = -1024
        NFramesToExclude = 3
        reuseCTMasked = false
        reuseCTRescaled = false
        reuseCarneyUmap = false 
        sessionDataCache
    end

    methods (Static)
        function assembleFdgAfterAC
            import mlsystem.* mlfourdfp.*;
            studyd = mlraichle.StudyData;            
            eSess = DirTool(fullfile(studyd.subjectsDir, 'HYGLY*'));
            for iSess = 1:length(eSess.fqdns)

                eVisit = DirTool(eSess.fqdns{iSess});
                for iVisit = 1:length(eVisit.fqdns)
                        
                    if (~isempty(regexp(eVisit.dns{iVisit}, '^V[1-2]$', 'match')))
                        fdgRawdata = sprintf('FDG_%s', eVisit.dns{iVisit});
                        pthAC = fullfile(eVisit.fqdns{iVisit}, [fdgRawdata '-AC'], '');
                        
                        if (isdir(pthAC))
                            rmdir(pthAC, 's');
                        end
                        
                        ensuredir(pthAC);
                        fprintf('MMRResolveBuilder.assembleFdgAfterAC:  working in -> %s\n', pthAC);                            
                        sessd = mlraichle.SessionData('studyData', studyd, ...
                                                      'sessionPath', eSess.fqdns{iSess}, ...
                                                      'tracer', 'FDG', ...
                                                      'vnumber', T4ResolveUtilities.visitNumber(eVisit.dns{iVisit}));
                        this = MMRResolveBuilder0('sessionData', sessd);  
                        firstFortranTimeFrame_ = this.firstFortranTimeFrame;                          
                        fdgACRevision = sessd.fdgACRevision('typ', 'fp');
                        fdgPrefix = sprintf('FDG_%s-LM-00-OP', eVisit.dns{iVisit});
                        fv = FourdfpVisitor;
                        eFrame = DirTool(fullfile(eVisit.fqdns{iVisit}, sprintf('%s-Converted-Frame*', fdgRawdata), ''));
                        for iFrame = 1:length(eFrame.fqdns)
                            try
                                pwd0 = pushd(eFrame.fqdns{iFrame});
                                fortranNumFrame = MMRResolveBuilder0.frameNumber(eFrame.dns{iFrame}, 1);
                                fdgFramename = MMRResolveBuilder0.frameFileprefix(fdgACRevision, fortranNumFrame);
                                fv.sif_4dfp(fdgPrefix);
                                fdgT4 = sprintf('%s_frame%i_to_resolved_t4', ...
                                                sessd.fdgNACRevision('typ', 'fp'), fortranNumFrame);
                                fqFdgT4 = fullfile(sessd.fdgT4Location, fdgT4);
                                fv.cropfrac_4dfp(0.5, fdgPrefix, fdgACRevision);
                                if (fortranNumFrame >= firstFortranTimeFrame_ && ...
                                    lexist(fqFdgT4, 'file'))
                                    fv.lns(fqFdgT4);
                                    fv.t4img_4dfp(fdgT4, fdgACRevision, 'options', ['-O' fdgACRevision]);                            
                                    fv.move_4dfp([fdgACRevision '_on_resolved'], ...                                
                                                 fullfile(pthAC, [fdgFramename '_on_resolved']));
                                else                           
                                    fv.move_4dfp(fdgACRevision, ...                                
                                                 fullfile(pthAC, [fdgFramename '_on_resolved']));
                                end
                                delete('*.4dfp.*')
                                delete([fdgACRevision '_frame*_to_resolved_t4']);
                                popd(pwd0);
                            catch ME
                                handwarning(ME);
                            end
                        end
                        pwd0 = pushd(fullfile(pthAC, ''));
                        ipr.dest = fdgACRevision;
                        ipr.frames = ones(1, length(eFrame.fqdns));
                        this.pasteFrames(ipr, 'on_resolved');
                        fv.imgblur_4dfp([fdgACRevision '_on_resolved'], 5.5);
                        delete(fullfile(pthAC, [fdgACRevision '_frame*_on_resolved.4dfp.*']));
                        popd(pwd0);
                    end
                end                
            end
        end
        function triggeringMoveLogs
            mlfourdfp.MMRResolveBuilder0.triggering('moveLogs___', ...
                'conditions', {'isVisit' 'isAnyTracer' 'isNAC' 'isNotEmpty'});
        end
        function moveConvertedToConvertedNAC
            import mlsystem.* mlfourdfp.*;
            studyd = mlraichle.StudyData;            
            eSess = DirTool(studyd.subjectsDir);
            for iSess = 1:length(eSess.fqdns)

                eVisit = DirTool(eSess.fqdns{iSess});
                for iVisit = 1:length(eVisit.fqdns)
                        
                    eTracer = DirTool(eVisit.fqdns{iVisit});
                    for iTracer = 1:length(eTracer.fqdns)

                        try
                            nacFold = [eTracer.fqdns{iTracer} '-NAC'];
                            if (~isempty(regexp(eTracer.dns{iTracer}, '\w+-Converted$', 'once')) && ~isdir(nacFold))
                                movefile(eTracer.fqdns{iTracer}, nacFold);
                            end
                        catch ME
                            handwarning(ME);
                        end
                    end
                end                
            end
        end
        function incomplete = scanForIncompleteE7
            pthPPG0 = pushd(mlraichle.RaichleRegistry.instance.subjectsDir);
            incomplete = {};
            
            import mlsystem.* mlfourdfp.*;
            dtSess = DirTool('HYGLY*');
            for eSess = 1:length(dtSess.dns)
                cd(dtSess.fqdns{eSess});
                dtVisit = DirTool('V*');
                for eVisit = 1:length(dtVisit.dns)
                    cd(dtVisit.fqdns{eVisit});
                    
                    dtTracer = DirTool(['*_' dtVisit.dns{eVisit}]);
                    for eTracer = 1:length(dtTracer.dns)
                        dnTracer = dtTracer.dns{eTracer};
                        fqdnTracer = dtTracer.fqdns{eTracer};
                        prefixTracer = strtok(dnTracer, '-'); 
                        if ( lstrfind(dnTracer, 'FDG') && ...
                            ~lexist(sprintf('%s-Converted-Frame63/%s-LM-00-OP_000_000.v', ...
                                            fqdnTracer, prefixTracer)))
                            incomplete = [incomplete [fqdnTracer '-Converted-Frame63']]; %#ok<AGROW>
                        end
                        if ( lstrfind(dnTracer, 'HO') && ...
                            ~MMRResolveBuilder0.isConverted(fqdnTracer))
                            incomplete = [incomplete [fqdnTracer '-Converted']]; %#ok<AGROW>
                        end
                        if ( lstrfind(dnTracer, 'OO') && ...
                            ~MMRResolveBuilder0.isConvertedAbs(fqdnTracer))
                            incomplete = [incomplete [fqdnTracer '-Converted-Abs']]; %#ok<AGROW>
                        end
                        if ( lstrfind(dnTracer, 'OC') && ...
                            ~lexist(sprintf('%s-Converted-NAC/%s-LM-00-OP_009_000.v', ...
                                            fqdnTracer, prefixTracer)))
                            incomplete = [incomplete [fqdnTracer '-Converted-NAC']]; %#ok<AGROW>
                        end
                    end
                    
                    dtTracer = DirTools('*-Converted', '*-Converted-Abs');
                    for eTracer = 1:length(dtTracer.dns)
                        dnTracer = dtTracer.dns{eTracer};
                        fqdnTracer = dtTracer.fqdns{eTracer};
                        prefixTracer = strtok(dtTracer.dns{eTracer}, '-');  
                        if ( lstrfind(dnTracer, 'HO') && ...
                            ~lexist(sprintf('%s/%s-LM-00/%s-LM-00-OP_057_000.v', ...
                                            fqdnTracer, prefixTracer, prefixTracer)))
                            incomplete = [incomplete fqdnTracer]; %#ok<AGROW>
                        end                        
                        if ( lstrfind(dnTracer, 'OO') && ...
                            ~lexist(sprintf('%s/%s-LM-00/%s-LM-00-OP_057_000.v', ...
                                            fqdnTracer, prefixTracer, prefixTracer)))
                            incomplete = [incomplete fqdnTracer]; %#ok<AGROW>
                        end
                    end
                end
            end
            popd(pthPPG0);
            
            incomplete = incomplete';
        end
        function tf = hasAbsInLog(fqfn)
            log = mlio.TextParser.load(fqfn);
            MARKER = 'command line: C:\Siemens\PET\bin.win64-VA20\e7_recon --abs';
            str = log.findFirstCell(MARKER);
            tf = lstrfind(str, '--abs');
        end
        function tf = isConverted(trpath)
            cpath  = sprintf('%s-Converted', trpath);
            tf1    = isdir(cpath);            
            lmpath = fullfile(cpath, [basename(trpath) '-LM-00'], '');
            try
                dt  = mlsystem.DirTool(fullfile(lmpath, 'log_e7_recon_*.txt'));
                tf2 = ~mlfourdfp.MMRResolveBuilder0.hasAbsInLog(dt.fqfns{end});
            catch
                tf2 = false;
            end
            
            tf = tf1 && tf2;
        end
        function tf = isConvertedAbs(trpath)
            cpath  = sprintf('%s-Converted-Abs', trpath);
            tf1    = isdir(cpath);            
            lmpath = fullfile(cpath, [basename(trpath) '-LM-00'], '');
            try
                dt  = mlsystem.DirTool(fullfile(lmpath, 'log_e7_recon_*.txt'));
                tf2 = mlfourdfp.MMRResolveBuilder0.hasAbsInLog(dt.fqfns{end});
            catch
                tf2 = false;
            end
            
            tf = tf1 && tf2;
        end
    end

	methods
 		function this = MMRResolveBuilder0(varargin)
 			%% MMRRESOLVEBUILDER0
 			%  Usage:  this = MMRResolveBuilder0()

 			this = this@mlfourdfp.T4ResolveBuilder0(varargin{:});            
        end        
        
        function [ctm,ic]          = buildCTMasked(this)
            %% CTMASKED
            %  @return ctm := this.sessionData.ctMasked('typ', 'fqfp')
            %  @return ic  := ctMasked as ImagingContext
            
            import mlfourd.*;
            mpr  = this.sessionData.mpr('typ', 'fqfp');
            ct   = this.sessionData.ct('typ', 'fqfp');
            ct_  = sprintf('%s_%s', ct, datestr(now, 30));
            ct__ = sprintf('%s%s', ct_, 'a');
            ctm  = this.sessionData.ctMasked('typ', 'fqfp');
            if (lexist(this.fourdfpImg(ctm)) && this.reuseCTMasked)
                ic = ImagingContext(ctm);
                return
            end
            
            [~,ctToMprT4] = this.CT2mpr_4dfp(ct);
            mprToCtT4 = this.buildVisitor.t4_inv(ctToMprT4);
            this.buildVisitor.imgblur_4dfp(mpr, 10);
            this.buildVisitor.maskimg_4dfp([this.sessionData.ct('typ', 'fp') '_on_mpr'], [mpr '_b100'], ct_, 'options', '-t5'); 
            this.buildVisitor.maskimg_4dfp(ct_, ct_, ct__, 'options', '-t50')
            this.buildVisitor.t4img_4dfp(mprToCtT4, ct__, 'out', ctm, 'options', ['-O' ct])
            delete([ct_ '.4dfp.*']);
            delete([ct__ '.4dfp.*']);
        end
        function imageOnSumt       = ctOnPetSumt(this, ct, petSumt)
            assert(lexist(this.fourdfpImg(ct)));
            assert(lexist(this.fourdfpImg(petSumt)));
            
            [this,ctToMprT4]  = this.CT2mpr_4dfp(ct);
            [this,petToMprT4] = this.petSumt2mpr(petSumt);
            mprToPetT4        = this.buildVisitor.t4_inv(petToMprT4);
            ctToPetT4         = this.buildVisitor.t4_mul(ctToMprT4, mprToPetT4);
            imageOnSumt       = this.buildVisitor.t4img_4dfp(ctToPetT4, ct, 'options', ['-O' petSumt]);
        end 
        function [this,ctToMprT4]  = CT2mpr_4dfp(this, ct)
            assert(lexist(this.fourdfpImg(ct), 'file'));
            ctToMprT4 = this.buildVisitor.filenameT4(ct, this.sessionData.mpr('typ', 'fp'));
            if (lexist(ctToMprT4, 'file'))
                return
            end
            
            this.buildVisitor.CT2mpr_4dfp(this.sessionData.mpr('typ', 'fqfp'), ct, 'options', ['-T' this.atlas('typ', 'fqfp')]);
        end 
        function                     ensureTracerSymlinks(this)
            sessd = this.sessionData;
            mprAtlT4 = [sessd.mpr('typ', 'fp') '_to_' sessd.atlas('typ', 'fp') '_t4'];
            fqMprAtlT4 = fullfile(sessd.mpr('typ', 'path'), mprAtlT4);
            
            assert(lexist(fqMprAtlT4, 'file'));
            assert(this.buildVisitor.lexist_4dfp(sessd.mpr('typ', 'fqfp')));
            assert(this.buildVisitor.lexist_4dfp(sessd.ct( 'typ', 'fqfp')));
            assert(isdir(sessd.tracerNACLocation));
            
            cd(sessd.tracerNACLocation);
            if (~lexist(mprAtlT4))
                this.buildVisitor.lns(fqMprAtlT4);
            end
            if (~lexist(sessd.mpr('typ', 'fn')))
                this.buildVisitor.lns_4dfp(sessd.mpr('typ', 'fqfp'));
            end
            if (~lexist(sessd.ct('typ', 'fn')))
                this.buildVisitor.lns_4dfp(sessd.ct('typ', 'fqfp'));
            end
        end
        function this              = loadSessionDataCache(this, varargin)
            ip = inputParser;
            addOptional(ip, 'tracers', {'OC' 'HO' 'OO'}, @iscell);
            parse(ip, varargin{:});
            
            cache = {};
            for s = 1:2
                for t = 1:length(ip.Results.tracers)
                    sessd = this.sessionData;
                    sessd.tracer = ip.Results.tracers{t};
                    sessd.snumber = s;
                    if (isdir(sessd.tracerListmodeLocation))
                        cache = [cache {sessd}]; %#ok<AGROW>
                    end
                end
            end
            this.sessionDataCache = cache;
        end
        function imageOnSumt       = mprOnPetSumt(this, mpr, petSumt)
            assert(lexist(this.fourdfpImg(mpr)));
            assert(lexist(this.fourdfpImg(petSumt)));
            
            [this,petToMprT4] = this.petSumt2mpr(petSumt);
            mprToPetT4        = this.buildVisitor.t4_inv(petToMprT4);
            imageOnSumt       = this.buildVisitor.t4img_4dfp(mprToPetT4, mpr, 'options', ['-O' petSumt]);
        end          
        function [this,petToMprT4] = petSumt2mpr(this, varargin)
            ip = inputParser;
            addRequired( ip, 'pet', @mlfourdfp.FourdfpVisitor.lexist_4dfp);
            addParameter(ip, 'pharynx', false, @islogical);
            parser(ip, varargin{:});
            pet = ip.Results.pet;
            
            mpr        = this.sessionData.mpr('typ', 'fqfp');
            petToMprT4 = [pet '_to_' mybasename(mpr) '_t4'];
            if (lexist(petToMprT4, 'file'))
                delete(petToMprT4);
            end
            
            this.buildVisitor.t4_inv( ...
                this.buildVisitor.sagittal_t4, 'out', petToMprT4);
            
            petToMprT4 = this.buildVisitor.align_multiSpectral( ...
                'dest', mpr, ...
                'source', pet, ...
                'destMask', 'none', ...
                'sourceMask', 'none', ...
                'destBlur', this.blurArg, ...
                'sourceBlur', this.blurArg, ...
                't40', petToMprT4, ...
                't4', petToMprT4);              
            mprToAtlT4 = [mybasename(mpr) '_to_' this.atlas('typ', 'fp') '_t4'];
            if (~lexist(mprToAtlT4, 'file'))     
                this = this.buildVisitor.msktgenMprage(mpr, this.atlas('typ', 'fp'));          
            end
            
            petToAtlT4 = [pet '_to_' this.atlas('typ', 'fp') '_t4'];
            this.buildVisitor.t4_mul(petToMprT4, mprToAtlT4, petToAtlT4);            
            this.buildVisitor.t4img_4dfp( ...
                petToAtlT4, pet, ...
                'out', [mybasename(pet), '_on_' this.atlas('typ', 'fp')], ...
                'options', ['-O' this.atlas('typ', 'fqfp')]);            
            this.buildVisitor.msktgen2_4dfp( ...
                pet, this.msktgenThresh, ...
                'options', ['-T' this.atlas('typ', 'fqfp')]);            
        end   
        function ct                = rescaleCT(this, varargin)
            ip = inputParser;
            addOptional(ip, 'ctMasked', this.sessionData.ctMasked('typ', 'fqfp'), @mlfourdfp.FourdfpVisitor.lexist_4dfp);
            parse(ip, varargin{:});      
            
            fqfn = this.sessionData.ctRescaled('typ', '4dfp.ifh');
            [~,fn] = fileparts(fqfn);  
            ct = myfileprefix(fqfn);
            if (lexist(fqfn, 'file') && lexist(fn, 'file') && this.reuseCTMasked && this.reuseCTRescaled)
                return
            end
            
            ic = mlfourd.ImagingContext([ip.Results.ctMasked '.4dfp.ifh']);
            ic = ic.numericalNiftid;
            ic = ic * this.ct_rescaleSlope + this.ct_rescaleIntercept;
            
            ic.noclobber = false;
            ic = ic.saveas(fqfn); %#ok<NASGU>
            this.buildVisitor.lns_4dfp(myfileprefix(fqfn));
        end  
        function this              = repairConvertedNAC(this, frame1st, frame2nd)
            %% REPAIRCONVERTEDNAC
            %  @param frame1st is numeric.
            %  @param frame2nd is numeric.            
            %  See also:  mlraichle.T4ResolveBuilder0.repairSingle
            
            sessd = this.sessionData;
            cd(sessd.tracerNAC('typ', 'path'));
            this.printv('repairConvertedNAC.pwd -> %s\n', pwd);
            this.ensureTracerSymlinks;
            this = this.repairSingle( ...
                frame1st, frame2nd, ...
                'dest', sprintf('%sv%ir%i', lower(sessd.tracer), sessd.vnumber, sessd.rnumber));
        end
        function this              = resolveConvertedNAC(this)
            %% RESOLVECONVERTEDNAC is the principle caller of resolve.
            
            sessd = this.sessionData;
            pwd0 = pushd(sessd.tracerNAC('typ', 'path'));
            this.printv('resolveConvertedNAC.pwd -> %s\n', pwd);
            this.ensureTracerSymlinks;
            this = this.resolve( ...
                'dest',      sessd.tracerNACRevision('typ', 'fp'), ... 
                'source',    sessd.tracerNAC('typ', 'fp'), ...
                'firstCrop', this.firstCrop, ...
                'frames',    this.frames);
            popd(pwd0);
        end
    end

    %% PROTECTED
    
    methods (Access = protected) 
        function fr   = firstFortranTimeFrame(this)
            NNativeFrames = this.readSize(this.sessionData.fdgNACRevision('typ', 'fqfp'));
            NUmapFrames   = this.readSize(this.sessionData.fdgNACResolved('typ', 'fqfp'));
            fr = NNativeFrames - NUmapFrames + 1;
        end
        function        moveLogs___(this, varargin)
            pwd0 = pushd(this.sessionData.tracerNACLocation);
            movefiles('*.log', this.logPath);
            movefiles('*.txt', this.logPath);
            popd(pwd0);
        end
    end
    
    %% PRIVATE
    
    methods (Static, Access = private)
        function this = triggering(varargin)
            
            studyd = mlraichle.StudyData;
            
            ip = inputParser;
            addRequired( ip, 'method', @ischar);
            addOptional( ip, 'args', {});
            addParameter(ip, 'subjectsDir', studyd.subjectsDir, @isdir);
            addParameter(ip, 'conditions', {}, @iscell);
            addParameter(ip, 'tag', '', @ischar);
            parse(ip, varargin{:});
            if (~strcmp(ip.Results.subjectsDir, studyd.subjectsDir))
                studyd.subjectsDir = ip.Results.subjectsDir;
            end
            
            import mlsystem.* mlfourdfp.* ;
            %%MMRResolveBuilder0.printv('triggering.ip.Results:  %s\n', struct2str(ip.Results));            
            eSess = DirTool(ip.Results.subjectsDir);
            %%MMRResolveBuilder0.printv('triggering.eSess:  %s\n', cell2str(eSess.dns));
            for iSess = 1:length(eSess.fqdns)

                eVisit = DirTool(eSess.fqdns{iSess});
                %%MMRResolveBuilder0.printv('triggering.eVisit:  %s\n', cell2str(eVisit.dns));
                for iVisit = 1:length(eVisit.fqdns)
                        
                    eTracer = DirTool(eVisit.fqdns{iVisit});
                    %%MMRResolveBuilder0.printv('triggering.eTracer:  %s\n', cell2str(eTracer.dns));
                    for iTracer = 1:length(eTracer.fqdns)

                        pth___ = eTracer.fqdns{iTracer};
                        %%MMRResolveBuilder0.printv('triggering.pth:  %s\n', pth___);
                        if (T4ResolveUtilities.matchesTag(pth___, ip.Results.tag) && ...
                            T4ResolveUtilities.pathConditions(pth___, ip.Results.conditions))
                            try
                                sessd = mlraichle.SessionData( ...
                                    'studyData',   studyd, ...
                                    'sessionPath', eSess.fqdns{iSess}, ...
                                    'snumber',     T4ResolveUtilities.scanNumber(eTracer.dns{iTracer}), ...
                                    'tracer',      T4ResolveUtilities.tracerPrefix(eTracer.dns{iTracer}), ...
                                    'vnumber',     T4ResolveUtilities.visitNumber(eVisit.dns{iVisit}));
                                this = MMRResolveBuilder0('sessionData', sessd);
                                this.(ip.Results.method)(ip.Results.args{:});   
                            catch ME
                                handwarning(ME);
                            end
                        end
                    end
                end                
            end            
        end
    end
    
    %  @deprecated
    
    methods (Hidden)
        function revertToLM00(this, nacPth)
            if (~isdir(nacPth))
                return
            end
            vPth = fileparts(nacPth);
            [~,vFold] = fileparts(vPth);
            tracerFold = [upper(this.sessionData.tracer) '_' vFold];
            lm00Pth = fullfile(vPth, [tracerFold '-Converted'], [tracerFold '-LM-00'], '');
            if (lexist(fullfile(nacPth, [lower(this.sessionData.tracer) lower(vFold) 'r2_resolved.4dfp.img'])))
                return
            end
            if (~isdir(lm00Pth))
                movefile(nacPth, lm00Pth);
            end
        end 
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

