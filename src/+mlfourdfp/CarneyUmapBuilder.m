classdef CarneyUmapBuilder < mlfourdfp.CTUmapBuilder
	%% CARNEYUMAPBUILDER  

	%  $Revision$
 	%  was created 05-Dec-2016 00:21:30
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.  Copyright 2017 John Joowon Lee.
 	   
    methods (Static)
        function buildUmapAll(varargin)
            ip = inputParser;
            addParameter(ip, 'tag', '', @ischar);
            parse(ip, varargin{:});

            import mlraichle.* mlsystem.*;
            setenv('PRINTV', '1');
            studyd = StudyData;            
            if (isempty(ip.Results.tag))
                tagString = 'HYGLY3*';
            else                
                tagString = [ip.Results.tag '*'];
            end
            
            eSess = DirTool(fullfile(studyd.subjectsDir, tagString));
            for iSess = 1:length(eSess.fqdns)

                try
                    pth = fullfile(eSess.fqdns{iSess}, '');
                    pwd0 = pushd(pth);
                    mlfourdfp.CarneyUmapBuilder.printv('buildUmapAll:  try pwd->%s\n', pwd);
                    sessd = SessionData( ...
                        'studyData',   studyd, ...
                        'sessionPath', eSess.fqdns{iSess}, ...
                        'ac',          false, ...
                        'tracer',      '');
                    if (mlfourdfp.FourdfpVisitor.lexist_4dfp(sessd.umapSynth('typ', 'fqfp')))
                        continue
                    end
                    this = mlfourdfp.CarneyUmapBuilder('sessionData', sessd);                        
                    this.keepForensics = true;
                    this.buildUmap;                        
                    popd(pwd0);
                catch ME
                    handwarning(ME);
                end
            end            
        end
    end
    
	methods 
        function [this,umap] = buildUmap(this, varargin)
            [this,umap] = this.buildCarneyUmap(varargin{:});
        end
        function [this,umap] = buildCarneyUmap(this, varargin)
            umap = [];
            if (this.isfinished)
                return
            end            
            pwd0 = pushd(this.sessionData.sessionPath);
            tracer0 = this.sessionData.tracer;
            this.sessionData_.tracer = 'FDG';
            this.convertUmapTo4dfp; % convert FDG_V*-Converted-NAC/FDG_V*-LM-00/FDG_V*-LM-00-umap.v
            this.sessionData_.tracer = tracer0;
            this.ensureSymlinkCTForms;
            ctm  = this.buildCTMasked; % ct_on_T1001 has excellent alignment
            ctm  = this.rescaleCT(ctm);
            umap = this.assembleCarneyUmap(ctm);
            umap = this.buildVisitor.imgblur_4dfp(umap, 4);
            %this.teardownBuildUmaps;
            popd(pwd0);
        end
        function [ctm,ic,ctToMprT4] = buildCTMasked(this)
            %% BUILDCTMASKED2 calls CT2mpr_4dfp
            %  @return ctm := this.sessionData.ctMasked('typ', 'fqfp')
            %  @return ic  := ctMasked as ImagingContext2 on MPR-space
            
            import mlfourd.*;
            mpr  = this.sessionData.mpr('typ', 'fqfp');
            ct   = this.sessionData.ct('typ', 'fqfp');
            ctm  = this.sessionData.ctMasked('typ', 'fqfp');
            ctd  = fullfile(this.sessionData.sessionPath, 'ct', '');
            
            assert(isdir(ctd));            
            pwd0 = pushd(fileparts(ctd));
            
            if (~lexist([ct '.4dfp.img'])) % disambiguate ctd and ct
                this.buildVisitor.dcm_to_4dfp(fullfile(ctd, '*.dcm'), 'base', 'ct', 'options', '-g');
            end
            
            [ctOnMpr,ctToMprT4] = this.CT2mpr_4dfp(ct, ...
                'log', sprintf('CarneyUmapBuilder_CT2mpr_4dfp_%s.log', mydatetimestr(now)));
            mprb = this.buildVisitor.imgblur_4dfp(mpr, 10);
            
            ct_  = sprintf('%s_%s', ct, mydatetimestr(now));
            this.buildVisitor.maskimg_4dfp(ctOnMpr, mprb, ct_, 'options', '-t5'); % in mpr-space
            this.buildVisitor.maskimg_4dfp(ct_, ct_, ctm, 'options', '-t50');
            ic = ImagingContext2(ctm);
            delete([ct_ '.4dfp.*']); % ct__ in mpr-space has best registration
            
            popd(pwd0);
        end
        function [this,umap] = buildPhantomUmap(this, varargin)
            pwd0 = pushd(this.sessionData.sessionPath);
            ip = inputParser;
            addOptional(ip, 'ctm', 'ctMasked', @lexist_4dfp);
            parse(ip, varargin{:});

            this.ensureSymlinkCTForms;
            ctr  = this.rescaleCT(ip.Results.ctm, 'ctOut', 'ctRescaled');
            umap = this.assembleCarneyUmap(ctr, 'umapSynth');
            umap = this.buildVisitor.imgblur_4dfp(umap, 4);
            popd(pwd0);
        end
        function umap        = testGroundTruth(this, ct)
            assert(lexist_4dfp(ct));
            ct = this.rescaleCT(ct, 'ctOut', [ct '_rescaled']);
            umap = this.assembleCarneyUmap(ct, [ct '_umap']);
            umap = this.buildVisitor.imgblur_4dfp(umap, 4);
        end
        function umaps       = umapsOpTracer(this)
            umaps = sprintf('umapsOp%sr%i', ...
                upperFirst(this.sessionData.tracer), this.sessionData.rnumber);
        end
        function fp          = umapsOpTracerFrame(this, fr)
            ipr = struct( ...
                'dest', this.umapsOpTracer, ...
                'currentIndex', fr);
            fp = this.fileprefixIndexed(ipr);
        end
		  
        %% E7 Utilities
        
        function this = convertUmapTo4dfp(this)
            pwd0 = pwd;
            if (~lexist(this.sessionData.tracerListmodeUmap('typ', '4dfp.img')))
                cd( this.sessionData.tracerListmodeUmap('typ', 'path'));
                this.buildVisitor.IFhdr_to_4dfp( ...
                    this.sessionData.tracerListmodeUmap('typ', 'v.hdr'));
                cd(pwd0);
            end
        end
        function this = convertUmapToE7Format(this, varargin)
            sessd = this.sessionData;
            sessd.rnumber = 1;
            this.mmrBuilder_ = mlsiemens.MMRBuilder('sessionData', this.sessionData);
            
            ip = inputParser;
            addOptional(ip, 'umap', ...
                fullfile(sessd.tracerNACLocation, ...
                    sprintf('%s_op_%s', sessd.umapSynth('typ', 'fp'), sessd.tracerNACRevision('typ', 'fp'))), ...
                @lexist_4dfp);
            addParameter(ip, 'zoom', this.mmrBuilder_.inverseCrop, @isnumeric);
            parse(ip, varargin{:});
            umap = ip.Results.umap;
            
            flipped = this.buildVisitor.flip_4dfp('z', umap);
            ic = mlfourd.ImagingContext2([flipped '.4dfp.hdr']);
            ic = ic.zoomed(ip.Results.zoom);
            ic.noclobber = false;
            ic.saveas([flipped '.4dfp.hdr']);
            movefile( ...
                sprintf('%s.4dfp.img', flipped), ...
                sprintf('%s.v',        umap), 'f');
            if (~this.keepForensics)
                delete(sprintf('%s.4dfp.*', flipped));
                delete(sprintf('%sfz.4dfp.*', umap));
                delete(sprintf('%s*.log', umap));
            end
            
            this.product_ = mlfourd.ImagingContext2(sprintf('%s.v', umap));
        end
        function this = convertUmapsToE7Format(this, umaps)
            assert(iscell(umaps));
            prodCell = {};
            for fr = 1:length(umaps)
                this = this.convertUmapToE7Format(umaps{fr});
                prodCell{fr} = this.product_.fqfilename; %#ok<AGROW>
            end
            
            this.product_ = prodCell;
        end  
        function this = repUmapToE7Format(this, umaps)
            assert(iscell(umaps));
            prodCell = {};
            for fr = 1:length(umaps)
                this = this.repUmapToE7Format__(umaps{fr});
                prodCell{fr} = this.product_.fqfilename; %#ok<AGROW>
            end
            
            this.product_ = prodCell;
        end 
        function this = repUmapToE7Format__(this, varargin)
            sessd = this.sessionData;
            sessd.rnumber = 1;
            this.mmrBuilder_ = mlsiemens.MMRBuilder('sessionData', this.sessionData);
            
            ip = inputParser;
            addOptional(ip, 'umap', ...
                fullfile(sessd.tracerNACLocation, ...
                    sprintf('%s_op_%s', sessd.umapSynth('typ', 'fp'), sessd.tracerNACRevision('typ', 'fp'))), ...
                @lexist_4dfp);            
            addParameter(ip, 'zoom', this.mmrBuilder_.inverseCrop, @isnumeric);
            parse(ip, varargin{:});
            umap = ip.Results.umap;
            
            flipped = this.buildVisitor.flip_4dfp('z', umap);
            ic = mlfourd.ImagingContext2([flipped '.4dfp.hdr']);
            ic = ic.zoomed(ip.Results.zoom);
            ic.noclobber = false;
            ic.saveas([flipped '.4dfp.hdr']);
            movefile( ...
                sprintf('%s.4dfp.img', flipped), ...
                sprintf('%s.v',        umap), 'f');
            if (~this.keepForensics)
                delete(sprintf('%s.4dfp.*', flipped));
                delete(sprintf('%sfz.4dfp.*', umap));
                delete(sprintf('%s*.log', umap));
            end
            
            this.product_ = mlfourd.ImagingContext2(sprintf('%s.v', umap));
        end 
        function        reconvertUmapsToE7Format(this)
            pwd0 = pushd(this.sessionData.fdgNACLocation);
            for fr = 1:length(this.indicesLogical)
                this.buildVisitor.extract_frame_4dfp(this.umapsOpTracer, fr);
                this.convertUmapToE7Format(this.umapsOpTracerFrame(fr));
                delete([this.umapsOpTracerFrame(fr) '.4dfp.*']);
                delete([this.umapsOpTracerFrame(fr) '_flipz.log']);
            end
            popd(pwd0);
        end
        
        %%
        
 		function this = CarneyUmapBuilder(varargin)
 			%% CARNEYUMAPBUILDER
 			%  Usage:  this = CarneyUmapBuilder()

 			this = this@mlfourdfp.CTUmapBuilder(varargin{:});
            this.finished_ = mlpipeline.Finished(this, ...
                'path', this.getLogPath, 'tag', lower(this.sessionData.tracer));
        end        
 	end 

    %% PROTECTED   
    
    properties (Access = protected)
        mmrBuilder_
    end
    
    methods (Access = protected)
        function umap = assembleCarneyUmap(this, varargin)
            %% ASSEMBLECARNEYUMAP follows Carney, et al. Med. Phys. 33(4) 2006 976-983.
            %  @param rescaledCT is the (fully-qualified) fileprefix of the rescaled CT.
            %  @param umap       is the (fully-qualified) fileprefix of the product umap.
            %  @returns umap product as f.-q. fileprefix.
            
            import mlfourdfp.*;
            ip = inputParser;
            addRequired(ip, 'rescaledCT', @FourdfpVisitor.lexist_4dfp);
            addOptional(ip, 'umap', this.sessionData.umapSynth('tracer', '', 'blurTag', '', 'typ', 'fqfp'), @ischar);
            parse(ip, varargin{:});            
            umap = ip.Results.umap;
            
            if (FourdfpVisitor.lexist_4dfp(umap) && this.reuseCarneyUmap)
                return
            end
            delete([umap '.4dfp.*']);
            delete([umap '_b40.4dfp.*']);
            delete([umap '.log']);
            ic = this.CarneyImagingContext(ip.Results.rescaledCT);
            ic.saveas([umap '.4dfp.hdr']);
        end
        function        ensureSymlinkCTForms(this)
            sd = this.sessionData;
            sp = sd.sessionPath;
            bv = this.buildVisitor;
            pwd0 = pushd(sd.sessionPath);
            if (~lexist_4dfp(sd.ct('typ', 'fp')))
                if (~lexist_4dfp(fullfile(sp, 'AC_CT')))
                    assert(lexist(fullfile(sp, 'AC_CT_series2.4dfp.hdr')))
                    bv.lns_4dfp(fullfile(sp, 'AC_CT_series2'), 'AC_CT');
                end
                bv.lns_4dfp(fullfile(sp, 'AC_CT'), sd.ct('typ', 'fp'));
            end
            popd(pwd0);
        end
        function umap = CarneyImagingContext(this, varargin)
            %% CARNEYIMAGINGCONTEXT follows Carney, et al. Med. Phys. 33(4) 2006 976-983.
            %  @param ct   is the (fully-qualified) fileprefix of the rescaled CT.
            %  @returns umap ImagingContext2.
            
            import mlfourdfp.*;
            ip = inputParser;
            addOptional(ip, 'ctRescaled', '', @FourdfpVisitor.lexist_4dfp);
            parse(ip, varargin{:});
            
            ct  = mlfourd.ImagingContext2([ip.Results.ctRescaled '.4dfp.hdr']);
            ct  = ct.selectNumericalTool;
            
            lowHU    = ct.uthresh(this.CarneyBP); 
            lowMask  = lowHU.binarized;
            highMask = lowMask.ones - lowMask;
            highHU   = ct .* highMask;
            
            lowHU = (lowHU + 1000) * 9.6e-5;
            lowHU =  lowHU .* lowMask;
            
            highHU = (highHU + 1000) * this.CarneyA + this.CarneyB;
            highHU =  highHU .* highMask;

            umap = lowHU + highHU;
            %umap = umap.blurred(4);
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
                    error('mlfourdfp:valueOutOfBounds', 'CarneyUmapBuilder.ct_kVp->%g is not supported', this.ct_kVp);
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
                    error('mlfourdfp:valueOutOfBounds', 'CarneyUmapBuilder.ct_kVp->%g is not supported', this.ct_kVp);
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
                    error('mlfourdfp:valueOutOfBounds', 'CarneyUmapBuilder.ct_kVp->%g is not supported', this.ct_kVp);
            end
        end 
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

