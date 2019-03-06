classdef CarneyUmapBuilder2 < mlfourdfp.AbstractUmapResolveBuilder
	%% CARNEYUMAPBUILDER2  
    %  N.B. CarneyImagingContext which requires flip.
    %  Refactoring:  pull up methods from AbstractUmapResolveBuilder and isolate from CompositeT4ResolveBuilder.

	%  $Revision$
 	%  was created 15-Nov-2018 15:49:50 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
        ct_kVp = 120
        reuseCarneyUmap = true
 	end

	methods 
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
            
            deleteExisting([umap '.4dfp.*']);
            deleteExisting([umap '_b*.4dfp.*']);
            deleteExisting([umap '.log']);
            ic = this.CarneyImagingContext(ip.Results.rescaledCT);
            ic.saveas([umap '.4dfp.hdr']);
        end   
        function fqfp = prepareBrainmaskMskt(this)
            fqfp = fullfile(this.sessionData.sessionPath, 'brainmask_mskt');
            if (~lexist_4dfp(fqfp))
                ic2 = mlfourd.ImagingContext2('brainmask.4dfp.hdr');
                ic2 = ic2.binarizeBlended(2*mlpet.Resources.instance.pointSpread);
                ic2 = ic2.selectNumericalTool;
                ic2 = ic2 * 1000;
                ic2.saveas([fqfp '.4dfp.hdr']);
            end
        end
        function        teardownBuildUmaps(this)
            this.teardownLogs;
            this.teardownT4s;
            deleteExisting(fullfile(this.sessionData.sessionPath, 'ctRescaled*.4dfp.*'));
            deleteExisting(fullfile(this.sessionData.sessionPath, 'ct_on_*.4dfp.*'));
            this.finished.markAsFinished( ...
                'path', this.logger.filepath, 'tag', [this.finished.tag '_' myclass(this) '_teardownBuildUmaps']); 
        end
		  
 		function this = CarneyUmapBuilder2(varargin)
 			this = this@mlfourdfp.AbstractUmapResolveBuilder(varargin{:});
            this.sessionData.attenuationCorrected = false;
            this.NRevisions = 2;
            this.finished_ = mlpipeline.Finished(this, ...
                'path', this.getLogPath, 'tag', lower(this.sessionData.tracer));
 		end
 	end 

    %% PROTECTED
    
    methods (Access = protected)
        function umap = CarneyImagingContext(this, varargin)
            %% CARNEYIMAGINGCONTEXT follows Carney, et al. Med. Phys. 33(4) 2006 976-983.
            %  @param ct   is the (fully-qualified) fileprefix of the rescaled CT.
            %  @returns umap ImagingContext2.
            
            import mlfourdfp.*;
            ip = inputParser;
            addOptional(ip, 'ctRescaled', '', @FourdfpVisitor.lexist_4dfp); % this.buildCTMasked
            parse(ip, varargin{:});
            
            ct  = mlfourd.ImagingContext2([ip.Results.ctRescaled '.4dfp.hdr']);
            %ct  = flip(ct, 1);
            
            lowHU    = ct.uthresh(this.CarneyBP); 
            lowMask  = lowHU.binarized;
            highMask = lowMask.ones - lowMask;
            highHU   = ct .* highMask;
            
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
                    error('mlfourdfp:valueOutOfBounds', 'CarneyUmapBuilder2.ct_kVp->%g is not supported', this.ct_kVp);
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
                    error('mlfourdfp:valueOutOfBounds', 'CarneyUmapBuilder2.ct_kVp->%g is not supported', this.ct_kVp);
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
                    error('mlfourdfp:valueOutOfBounds', 'CarneyUmapBuilder2.ct_kVp->%g is not supported', this.ct_kVp);
            end
        end 
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

