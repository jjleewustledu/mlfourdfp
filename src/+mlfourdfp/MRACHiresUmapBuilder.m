classdef MRACHiresUmapBuilder < mlfourdfp.MRUmapBuilder
	%% MRACHIRESUMAPBUILDER 
    
    %  [jjlee@pascal umaps] history
    %  ...
    %  1610  Apr 25 00:02:40: dcm2niix -f test_dcm2niix -o `pwd` -z y Head_MRAC_Brain_HiRes_in_UMAP_DT20181213121340.482500
    %  1611  Apr 25 00:02:44: ls
    %  1612  Apr 25 00:14:27: pwd
    %  1613  Apr 25 00:14:31: ls ..
    %  1614  Apr 25 00:14:40: ls ../mri
    %  1615  Apr 25 00:15:14: mri_convert ../mri/T1.mgz T1.nii
    %  1616  Apr 25 00:15:26: nifti_4dfp T1.nii T1.4dfp.hdr
    %  1617  Apr 25 00:15:32: nifti_4dfp -4 T1.nii T1.4dfp.hdr
    %  1618  Apr 25 00:15:35: ls
    %  1619  Apr 25 00:15:58: gunzip test_dcm2niix_e2_ph.nii.gz
    %  1620  Apr 25 00:16:14: nifti_4dfp -4 test_dcm2niix_e2_ph.nii test_dcm2niix_e2_ph.4dfp.hdr
    %  1621  Apr 25 00:16:16: ls
    %  1622  Apr 25 00:16:43: CT2mpr_4dfp
    %  1623  Apr 25 00:18:52: emacs -nw 
    %  1624  Apr 25 00:19:51: mpr2atl_4dfp
    %  1625  Apr 25 00:19:57: mpr2atl1_4dfp
    %  1626  Apr 25 00:21:16: mpr2atl1_4dfp T1 -T$REFDIR/TRIO_Y_NDC
    %  1627  Apr 25 00:21:34: ls
    %  1628  Apr 25 00:21:52: CT2mpr_4dfp T1 test_dcm2niix_e2_ph -T$REFDIR/TRIO_Y_NDC -m    
    %  ...

	%  $Revision$
 	%  was created 23-Apr-2019 12:01:22 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.5.0.1067069 (R2018b) Update 4 for MACI64.  Copyright 2019 John Joowon Lee.
 	
	properties (Constant)
 		PREFIX = 'UMAP'
        REUSE_UMAP = true
    end
    
    methods (Static)
        function out = flirt(varargin)
            ip = inputParser;
            addParameter(ip, 'in', [], @ischar)
            addParameter(ip, 'ref', [], @ischar)
            addParameter(ip, 'out', '', @ischar)
            addParameter(ip, 'omat', '', @ischar)
            addParameter(ip, 'bins', 255, @isnumeric)
            addParameter(ip, 'cost', 'corratio', @ischar)
            addParameter(ip, 'dof', 12, @isnumeric)
            addParameter(ip, 'interp', 'trilinear', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            assert(isfile([ipr.in '.nii']) || isfile([ipr.in '.nii.gz']))
            if isempty(ipr.out)
                ipr.out = [mybasename(ipr.in) '_on_' mybasename(ipr.ref)];
            end
            if isempty(ipr.omat)
                ipr.omat = [mybasename(ipr.in) '_on_' mybasename(ipr.ref) '.mat'];
            end
            system(sprintf(['/usr/local/fsl/bin/flirt -in %s -ref %s -out %s -omat %s ' ...
                            '-bins %i -cost %s -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof %i  -interp %s'], ...
                ipr.in, ipr.ref, ipr.out, ipr.omat, ipr.bins, ipr.cost, ipr.dof, ipr.interp))            
            if isfile([ipr.out '.nii.gz'])
                gunzip([ipr.out '.nii.gz'])
            end   
            out = ipr.out;
        end
    end

	methods
        function out = buildUmap(this)
            umap = this.dcm_to_UMAP; % e.g., sessionPath/umaps/UMAP_DT20180822112042.600000.nii
            sesp = this.sessionData.sessionPath;
            out = this.flirt( ...
                'in', umap, ...
                'ref', fullfile(sesp, 'T1'), ...
                'out', fullfile(sesp, this.umapSynthOpT1001('blurTag', '', 'typ', 'fp')));
            out = mlfourd.ImagingContext2(myniftiprefix(out));
            out.selectFourdfpTool();
            out.save();
        end        
        function fqfp = dcm_to_UMAP(this)
            %% transforms MR AC as DICOM to UMAP as 4dfp; works in this.umapPath
            %  @return fullfile(this.umapPath, this.fileprefix)
            
            fqfp = fullfile(this.umapPath, this.fileprefix());
            if this.REUSE_UMAP && isfile([fqfp '.4dfp.hdr'])
                return
            end
                
            pwd0 = pushd(this.umapPath);
            system(sprintf('dcm2niix -f %s -o %s -z n %s', ...
                this.filename_nii_gz, ...
                this.umapPath, ...
                mybasename(this.umapDicomPath)))
            ic2 = mlfourd.ImagingContext(this.filename_nii_gz);
            ic2 = ic2.selectNumericalTool;
            ic2 = ic2 / 1000;
            ic2.save            
            system(sprintf('nifti_4dfp -4 %s %s', gunzip(this.filename_nii_gz), this.filename_4dfp_hdr));
            popd(pwd0);
        end
        function [umapOnMpr,umapToMprT4] = CT2mpr_4dfp(this, umap, varargin)
            %% builds MPR and Atlas T4s de novo.
            
            assert(lexist(this.fourdfpImg(umap), 'file'), ...
                'mlfourdfp:RuntimeError', ...
                'MRACHiresUmapBuilder.CT2mpr_4dfp could not find %s', this.fourdfpImg(umap));
            mpr = this.sessionData.mpr('typ', 'fqfp');
            pth = fileparts(mpr);
            umapToMprT4 = fullfile(pth, this.buildVisitor.filenameT4(mybasename(umap), mybasename(mpr))); 
            umapOnMpr   = fullfile(pth, [mybasename(umap) '_on_' mybasename(mpr)]);   
            if (~lexist([mpr '_to_' this.atlas('typ','fp') '_t4']))
                this.buildVisitor.mpr2atl_4dfp(mpr);
            end
            if (~lexist(this.fourdfpImg(umapOnMpr)))
                umapOnMpr = this.buildVisitor.CT2mpr_4dfp(mpr, umap, ...
                    'options', ['-T' this.atlas('typ','fqfp')], varargin{:});
            end
            assert(lexist(umapToMprT4, 'file'));        
        end
        function s = mrSeriesLabel(~)
            s = 'Head_MRAC_Brain_HiRes_in_UMAP';
        end
        function teardownBuildUmaps(this)
            this.teardownLogs;
            this.teardownT4s;
            deleteExisting(fullfile(this.sessionData.sessionPath, [this.fileprefix '.4dfp.*']));
            deleteExisting(fullfile(this.sessionData.sessionPath, [this.fileprefix '_on_*.4dfp.*']));
            this.finished.markAsFinished( ...
                'path', this.logger.filepath, 'tag', [this.finished.tag '_' myclass(this) '_teardownBuildUmaps']); 
        end       
        
 		function this = MRACHiresUmapBuilder(varargin)
 			%% MRACHIRESUMAPBUILDER
 			%  @param .

 			this = this@mlfourdfp.MRUmapBuilder(varargin{:});
 		end
 	end 
    
    %% PROTECTED
    
    methods (Access = protected)
        function f  = filename_nii_gz(this)
            f = sprintf('%s_e2_ph.nii.gz', this.fileprefix());
        end
        function f  = filename_4dfp_hdr(this)
            f = [this.fileprefix() '.4dfp.hdr'];
        end
        function f  = fileprefix(this)
            f = sprintf('%s_%s', this.PREFIX, this.DTstring());
        end
        function dt = DTstring(this, varargin)
            ip = inputParser;
            addOptional(ip, 's', mybasename(this.umapDicomPath), @ischar);
            parse(ip);
            r = regexp(ip.Results.s, '(DT\d{14}.\d{6})', 'match');
            assert(~isempty(r));
            dt = r{1};
        end 
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

