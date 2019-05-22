classdef (Abstract) AbstractBuilder < mlpipeline.AbstractBuilder
	%% ABSTRACTBUILDER  

	%  $Revision$
 	%  was created 14-Nov-2018 15:01:10 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
 		
 	end

    methods (Static)
        function fn    = fourdfpHdr(fp)
            fn = [myfileprefix(fp) '.4dfp.hdr'];
        end 
        function fn    = fourdfpIfh(fp)
            fn = [myfileprefix(fp) '.4dfp.ifh'];
        end 
        function fn    = fourdfpImg(fp)
            fn = [myfileprefix(fp) '.4dfp.img'];
        end 
        function fn    = fourdfpImgRec(fp)
            fn = [myfileprefix(fp) '.4dfp.img.rec'];
        end
        function fn    = fslchfiletype(fn, varargin)
            ip = inputParser;
            addRequired(ip, 'fn', @(x) lexist(x, 'file'));
            addOptional(ip, 'type', 'NIFTI_GZ', @ischar);
            parse(ip, fn, varargin{:});
            
            fprintf('mlfourdfp.AbstractBuilder.fslchfiletype is working on %s\n', ip.Results.fn);
            mlbash(sprintf('fslchfiletype %s %s', 'NIFTI_GZ', ip.Results.fn));
            [p,f] = myfileparts(fn);
            fn = fullfile(p, [f mlfourd.NIfTIInfo.FILETYPE_EXT]);
        end
        function fn    = mri_convert(fn, varargin)
            %% MRI_CONVERT
            %  @param fn is the source possessing a filename extension recognized by mri_convert
            %  @param fn is the destination, also recognized by mri_convert.  Optional.  Default is [fileprefix(fn) '.nii.gz'] 
            
            import mlfourdfp.AbstractBuilder.niigzFilename;
            ip = inputParser;
            addRequired(ip, 'fn',  @(x) lexist(x, 'file'));
            addOptional(ip, 'fn2', niigzFilename(fn), @ischar);
            parse(ip, fn, varargin{:});            
            
            fprintf('mlfourdfp.AbstractBuilder.mri_convert is working on %s\n', ip.Results.fn);
            mlbash(sprintf('mri_convert %s %s', ip.Results.fn, ip.Results.fn2));
            fn = ip.Results.fn2;
        end
        function [s,r] = nifti_4dfp_4(varargin)
            vtor = mlfourdfp.FourdfpVisitor;
            [s,r] = vtor.nifti_4dfp_4(varargin{:});
        end
        function [s,r] = nifti_4dfp_n(varargin)
            vtor = mlfourdfp.FourdfpVisitor;
            [s,r] = vtor.nifti_4dfp_n(varargin{:});
        end
        function fn    = niiFilename(fn)
            [p,f] = myfileparts(fn);
            fn = fullfile(p, [f '.nii']);
        end
        function fn    = niigzFilename(fn)
            [p,f] = myfileparts(fn);
            fn = fullfile(p, [f '.nii.gz']);
        end
    end
    
	methods 
        function [s,r] = ensure4dfp(this, varargin)
            %% ENSURE4DFP
            %  @param filename is any string descriptor found in an existing file on the filesystem;
            %  ensureNifti will search for files with extensions .4dfp.hdr.
            %  @returns s, the bash status; r, any bash messages.  ensure4dfp ensures files are .4dfp.hdr.            
            
            ip = inputParser;
            addRequired(ip, 'filename', @ischar);
            parse(ip, varargin{:});
            
            s = 0; r = '';
            if (2 == exist(ip.Results.filename, 'file'))
                if (lstrfind(ip.Results.filename, '.4dfp'))
                    return
                end
                if (lstrfind(ip.Results.filename, '.mgz'))
                    fp = myfileprefix(ip.Results.filename);
                    [s,r] = mlbash(sprintf('mri_convert %s.mgz %s.nii', fp, fp)); %#ok<ASGLU>
                    [s,r] = this.ensure4dfp([fp '.nii']);
                    return
                end
                [s,r] = this.buildVisitor.nifti_4dfp_4(myfileprefix(ip.Results.filename));
                assert(lexist(myfilename(ip.Results.filename, '.4dfp.hdr'), 'file'));
                return
            end
            if (2 == exist([ip.Results.filename '.4dfp.hdr'], 'file'))
                return
            end
            if (2 == exist([ip.Results.filename '.nii'], 'file'))
                [s,r] = this.ensure4dfp([ip.Results.filename '.nii']);
                return
            end
            if (2 == exist([ip.Results.filename '.nii.gz'], 'file'))
                [s,r] = this.ensure4dfp([ip.Results.filename '.nii.gz']);
                return
            end      
            if (2 == exist([ip.Results.filename '.mgz'], 'file'))
                [s,r] = mlbash(sprintf('mri_convert %s.mgz %s.nii', ip.Results.filename, ip.Results.filename)); %#ok<ASGLU>
                [s,r] = this.ensure4dfp([ip.Results.filename '.nii']);
                return
            end            
            error('mlfourdfp:fileNotFound', ...
                  'T4ResolveBuilder.ensureNifti could not find files of form %s', ip.Results.filename);     
        end
        function [s,r] = ensureNifti(this, varargin)
            %% ENSURENIFTI
            %  @param filename is any string descriptor found in an existing file on the filesystem;
            %  ensureNifti will search for files with extensions .nii, .nii.gz or .4dfp.hdr.
            %  @returns s, the bash status; r, any bash messages.  ensureNifti ensures files are .nii.gz.
            
            ip = inputParser;
            addRequired(ip, 'filename', @ischar);
            parse(ip, varargin{:});
            
            s = 0; r = '';
            if (2 == exist(ip.Results.filename, 'file'))
                if (lstrfind(ip.Results.filename, '.nii'))
                    return
                end
                if (lstrfind(ip.Results.filename, '.mgz'))
                    fp = myfileprefix(ip.Results.filename);
                    [s,r] = mlbash(sprintf('mri_convert %s.mgz %s.nii.gz', fp, fp));
                    return
                end
                [s,r] = this.buildVisitor.nifti_4dfp_n(myfileprefix(ip.Results.filename));
                assert(lexist(myfilename(ip.Results.filename), 'file'));
                return
            end
            if (2 == exist([ip.Results.filename '.nii'], 'file'))
                [s,r] = this.ensureNifti([ip.Results.filename '.nii']);
                return
            end
            if (2 == exist([ip.Results.filename '.nii.gz'], 'file'))
                return
            end
            if (2 == exist([ip.Results.filename '.mgz'], 'file'))
                [s,r] = mlbash(sprintf('mri_convert %s.mgz %s.nii.gz', ip.Results.filename, ip.Results.filename));
                return
            end    
            if (2 == exist([ip.Results.filename '.4dfp.hdr'], 'file'))
                if (2 == exist([ip.Results.filename '.nii'], 'file'))
                    return
                end
                [s,r] = this.ensureNifti([ip.Results.filename '.4dfp.hdr']);
                return
            end
            error('mlfourdfp:fileNotFound', ...
                  'T4ResolveBuilder.ensureNifti could not find files of form %s', ip.Results.filename);            
        end
        function fps   = ensureSafeFileprefix(this, varargin)
            fps = this.buildVisitor.ensureSafeFileprefix(varargin{:});
        end        
		  
 		function this = AbstractBuilder(varargin)
 			%% ABSTRACTBUILDER
 			%  @param .

 			this = this@mlpipeline.AbstractBuilder(varargin{:});
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

