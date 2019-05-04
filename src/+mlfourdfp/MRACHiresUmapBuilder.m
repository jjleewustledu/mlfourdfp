classdef MRACHiresUmapBuilder < mlfourdfp.MRUmapBuider
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
 	
	properties
 		
 	end

	methods
        function umap = buildUmap(this)
            pct = this.dcm_to_pseudoct; 
            pct = this.rescale_pseudoct(pct);
            pct_on_mpr = this.CT2mpr_4dfp(this, pct);
            carneyBuilder = mlfourdfp.CarneyUmapBuilder2('sessionData', this.sessionData);
            umap = carneyBuilder.assembleCarneyUmap(pct_on_mpr);
        end        
        
        function fqfn = CT2mpr_4dfp(this, pct)
            ctBuilder = mlfourdfp.PseudoCTBuilder('sessionData', this.sessionData);
            fqfp = ctBuilder.CT2mpr_4dfp(pct);
            fqfn = this.fourdfpImg(fqfp);
        end
        function fqfn = dcm_to_pseudoct(this)
            %% transforms MR AC as DICOM to pseudo-CT as 4dfp
            
            pwd0 = pushd(this.umapPath);
            system(sprintf('dcm2niix -f %s -o %s -z y %s', ...
                this.filename_nii_gz, ...
                this.umapPath, ...
                mybasename(this.umapDicomPath)))
            system(sprintf('nifti_4dfp -4 %s %s', gunzip(this.filename_nii_gz), this.filename_4dfp_hdr));
            fqfn = fullfile(this.umapPath, this.filename_4dfp_hdr);
            popd(pwd0);
        end
        function s    = mrSeriesLabel(~)
            s = 'Head_MRAC_Brain_HiRes_in_UMAP';
        end
        function fqfn = rescale_pseudoct(pct)
            ctBuilder = mlfourdfp.PseudoCTBuilder('sessionData', this.sessionData);
            fqfp = [myfileprefix(pct) '_rescaled'];
            ctBuilder.rescaleCT(myfileprefix(pct), 'ctOut', fqfp)
            fqfn = this.fourdfpImg(fqfp);
        end
        function        teardownBuildUmaps(this)
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

 			this = this@mlfourdfp.MRUmapBuider(varargin{:});
 		end
 	end 
    
    %% PROTECTED
    
    methods (Access = protected)
        function f  = filename_nii_gz(this)
            f = sprintf('pseudoct_%s_e2_ph.nii.gz', this.DTstring());
        end
        function f  = filename_4dfp_hdr(this)
            f = sprintf('pseudoct_%s.4dfp.hdr', this.DTstring());
        end
        function f  = fileprefix(this)
            f = sprintf('pseudoct_%s', this.DTstring());
        end
        function dt = DTstring(~, varargin)
            ip = inputParser;
            addOptional(ip, 's', mybasename(this.umapDicomPath), @ischar);
            parse(ip);
            r = regexp(ip.Results.s, 'Head_MRAC_Brain_HiRes_in_UMAP_(DT\d{14}.\d{6})', 'match');
            assert(~isempty(r));
            dt = r{1};
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

