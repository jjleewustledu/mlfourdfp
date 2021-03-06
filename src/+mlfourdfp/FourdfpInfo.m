classdef FourdfpInfo < mlfourd.Analyze75Info
	%% FOURDFPINFO 
    %  See also ${TRANSFER}/4dfp_format.txt
    
    %% 4dfp OVERVIEW
    %  -------------
    % The 4dfp (4-dimensional floating point) format was designed for functional neuroimaging.
    % The four dimensions typically correspond to x, y, z, and time.
    % Three dimensional structural images can be represented in 4dfp format by setting the
    % depth of the fourth dimension to 1.
    % 
    % The voxel data are stored in the form of a binary image as one UNIX file.
    % Consequently, 4dfp images may be directly loaded and viewed using IDL, matlab, anaylue-avw, etc.
    % Information critical to interpreting the binary data (e.g., orientation,
    % image dimensions, voxel dimensions) are stored in separate header file[s].
    % The 4dfp UNIX file name convention is demonstrated below.
    % 
    % --------------------------------------------------------------
    % <filename>.4dfp.img		# binary float voxel data
    % <filename>.4dfp.ifh		# interfile header (ASCII text)
    % <filename>.4dfp.hdr		# ANALYZE 7.5 header (binary)
    % <filename>.4dfp.img.rec		# creation history
    % --------------------------------------------------------------
    % Please note that <filename> is any valid filename string.
    % 
    % All 4dfp based image analysis programs used at the Washington University School of Medicine
    % Neuroimaging Laboratory (NIL) read/write interfile headers.
    % The minimal 4dfp format is comprised of <filename>.4dfp.img + <filename>.4dfp.ifh.
    % To view 4dfp images using analyze_avw it is necessary to create an
    % ANALYZE header file, i.e., <filename>.4dfp.hdr. A utility for this purpose (ifh2hdr) is available.
    % All NIL image analysis programs maintain an additional file,
    % <filename>.4dfp.img.rec, which records the sequence of steps by which the image was created.
     
    %% 4dfp image data
    %  ---------------
    % The 4dfp file format was developed to manage human head images acquired with a Siemens MRI scanner.
    % The imato4dfp utility converts Siemens slice-based image (*.ima) files into 4dfp format.
    % This is accomplished by extracting (without reordering) the stored pixel values, converting short int to
    % to float and writing the results to the 4dfp image file.
    % Within each 3D volume the Siemens *.ima files are read in order of decreasing image number.
    % 
    % Certain 4dfp conventions, especially regarding image orientation, reflect the
    % Siemens Numaris operating system. For all acquisitions, including
    % oblique and double oblique, Numaris determines the principal orientation.
    % It is this orientation which appears in the interfile header and to which following table refers.
    % The orientation-dependent axis flips
    % required to display 4dfp image data in conventional radiologic orientation are tabulated below.
    % It is assumed that the display is written left-to-right and bottom-to-top as in analyze_avw.
    % It is also assumed that the image was acquired with the
    % subject positioned in the scanner head-first and supine.
    % 
    % --------------------------------------------------------------------
    % orientation		flip axes	first voxel (after flipping)
    % --------------------------------------------------------------------
    % 2 (axial)		y		right occipital skull-base
    % 3 (coronal)	 	y z		right skull-base occipital
    % 4 (sagittal)		x y z		occipital skull-base right
    % --------------------------------------------------------------------
    % 
    % The above enumerated axis flips may be effected in analyze_avw at 4dfp image load time.
    % Alternatively, 4dfp format may be converted to ANALYZE 7.5 format using the
    % 4dfptoanalyze utility which automatically performs the indicated axis flips as it
    % converts the voxel value number format from float to short int.
    % 4dfp image data loaded into analyze_avw as prescribed above will be consistently
    % resliced by analyze_avw. That is, all Siemens acquisition orientations will be
    % consistently displayed in all analyze orientations (transverse, coronal, sagittal).
    % 
    % Please note that x, y, and z denote voxel indices as ordered in memory.
    % That is, if the x, y, and z indices run from 0 to nx-1,
    % 0 to ny-1, and 0 to nz-1, then coordinate (x,y,z) is stored relative to the
    % first voxel of each frame at offset x + nx*(y + ny*z). 
     
    %% interfile header
    %  ----------------
    % The following is a listing of a 4dfp interfile header file (vm6c_b1.4dfp.ifh).
    % The image data (vm6c_b1.4dfp.img) were acquired in one 128 frame fMRI run.
    % Each frame has dimensions 64 x 64 x 18, The acquired voxels are 3 mm cubic.
    % 
    % ----------------------------------------------------------------------------------------
    % INTERFILE 			:=
    % version of keys			:= 3.3
    % image modality			:= mri
    % originating system		:= MAGNETOM VISION
    % conversion program		:= imato4dfp2
    % program version			:= $Id: imato4dfp2.c,v 1.12 2000/05/05 00:56:18 avi Exp $
    % original institution		:= Washington University
    % name of data file		:= vm6c_b1
    % patient ID 			:= vm6c
    % date				:= 22-FEB-1999
    % number format			:= float
    % number of bytes per pixel	:= 4
    % orientation			:= 3
    % time series flag		:= 1
    % number of dimensions		:= 4
    % scaling factor (mm/pixel) [1]	:= 3.000000
    % scaling factor (mm/pixel) [2]	:= 3.000000
    % scaling factor (mm/pixel) [3]	:= 3.000000
    % slice thickness (mm/pixel)	:= 3.000000
    % matrix size [1]			:= 64
    % matrix size [2]			:= 64
    % matrix size [3]			:= 18
    % matrix size [4]			:= 128
    % global minimum			:= 0
    % global maximum			:= 4095
    % mri parameter file name		:= Initialized by sequence
    % mri sequence file name		:= /usr/users/tec/nbea_uc_tg2.ekc
    % mri sequence description	:= ep_fid
    % -----------------------------------------------------------------------------------------
    % 
    % Various image analysis programs may make use of additional interfile header fields.
    % The minimal set of fields required to interpret the voxel data is listed below.
    % -------------------------------------------------
    % number format			:= float
    % number of bytes per pixel	:= 4
    % orientation			:= 3
    % number of dimensions		:= 4
    % scaling factor (mm/pixel) [1]	:= 3.000000
    % scaling factor (mm/pixel) [2]	:= 3.000000
    % scaling factor (mm/pixel) [3]	:= 3.000000
    % matrix size [1]			:= 64
    % matrix size [2]			:= 64
    % matrix size [3]			:= 18
    % matrix size [4]			:= 128
    % -------------------------------------------------
     
    %% rec file
    %  --------
    % The rec file format was designed to capture the creation history of each
    % particular 4dfp image. This is accomplished automatically provided that each
    % UNIX executable which creates 4dfp output also produces a corresponding rec file.
    % Rec files are ASCII text with the following format.
    % 
    % --------------------------------------------------------
    % rec <filename>.4dfp.img `date` `user`
    % UNIX command line which created <filename>.4dfp.img
    % rcs $Id$ (program revision code)
    % image/program specific processing information
    % ...
    % rec file[s] corresponding to antecedent input 4dfp images
    % endrec `date` `user`
    % --------------------------------------------------------
    % 
    % The critical feature of the rec file convention is
    % inclusion of antecedent rec files at all stages of processing. 
    % It follows that rec files corresponding to averaged images may grow large.
    % The key words "rec" (first field of first line) and "endrec"
    % (first field of last line) guarantee secure parsing of the accumulated
    % processing history. The brec (beautify rec file) utility parses rec files
    % and writes to stdout a more easily readable version of the text.
    % The following is a listing of the rec file corresponding to the above illustrated
    % interfile header.
    % 
    % ---------------------------------------------------------------------------------------------
    % rec vm6c_b1.4dfp.img  Thu May 18 17:15:18 2000  avi
    % /data/petsun4/data1/solaris/imato4dfp2 -fy /data/petsun23/vm6c/siem_im/bold1/5250 7 7 vm6c_b1 
    % $Id: imato4dfp2.c,v 1.12 2000/05/05 00:56:18 avi Exp $
    % patient_id:		vm6c
    % institution:		Washington University
    % manufacturer_model:	MAGNETOM VISION
    % parameter_file_name:	Initialized by sequence
    % sequence_file_name:	/usr/users/tec/nbea_uc_tg2.ekc
    % sequence_description:	ep_fid   90	TR    135.2	TE   37.0/1
    % tilts:			Cor>Tra -12	           
    % 4dfp_dimensions:	64        64        18        128       
    % voxel_dimensions:	3.000000  3.000000  3.000000  
    % scan_date:		22-FEB-1999
    % scan_time:		14:06:33-14:06:33
    % endrec Thu May 18 17:15:18 2000  avi
    % ---------------------------------------------------------------------------------------------
    % 
    % The following is vm6c_b1.4dfp.img.rec filtered through brec.
    % ----------------------------------------------------------------------------------------------------
    % 1rec vm6c_b1.4dfp.img  Thu May 18 17:15:18 2000  avi
    % 1      /data/petsun4/data1/solaris/imato4dfp2 -fy /data/petsun23/vm6c/siem_im/bold1/5250 7 7 vm6c_b1 
    % 1      $Id: imato4dfp2.c,v 1.12 2000/05/05 00:56:18 avi Exp $
    % 1      patient_id:             vm6c
    % 1      institution:            Washington University
    % 1      manufacturer_model:     MAGNETOM VISION
    % 1      parameter_file_name:    Initialized by sequence
    % 1      sequence_file_name:     /usr/users/tec/nbea_uc_tg2.ekc
    % 1      sequence_description:   ep_fid   90     TR    135.2     TE   37.0/1
    % 1      tilts:                  Cor>Tra -12                
    % 1      4dfp_dimensions:        64        64        18        128       
    % 1      voxel_dimensions:       3.000000  3.000000  3.000000  
    % 1      scan_date:              22-FEB-1999
    % 1      scan_time:              14:06:33-14:06:33
    % 1endrec Thu May 18 17:15:18 2000  avi
    % ----------------------------------------------------------------------------------------------------
    % 
    % The following is the rec file corresponding to the same 4dfp image data after
    % passage throgh the programs rmspike_4dfp and deband_4dfp (spike artifact removal and
    % even/odd slice banding compensation).
    % ----------------------------------------------------------------------------------------------------
    % rec vm6c_b1_rmsp_dbnd.4dfp.img  Thu May 18 17:16:23 2000  avi
    % /data/petsun4/data1/solaris/deband_4dfp -n4 vm6c_b1_rmsp 
    % $Id: deband_4dfp.c,v 1.8 1999/11/20 00:55:49 avi Exp $
    % Frame          1 slice multipliers: even=0.837060 odd=1.162940
    % Frame          2 slice multipliers: even=0.997099 odd=1.002901
    % Frame          3 slice multipliers: even=0.985484 odd=1.014516
    % Frame          4 slice multipliers: even=0.986583 odd=1.013417
    % Functional frame slice multipliers: even=0.986982 odd=1.013018
    % rec vm6c_b1_rmsp.4dfp.img  Thu May 18 17:16:13 2000 avi
    % /data/petsun4/data1/solaris/rmspike_4dfp -n4 -x33 vm6c_b1 
    % $Header: /data/petsun4/src_solaris/rmspike_4dfp/RCS/rmspike_4dfp.c,v 2.6 1997/05/23 00:49:24 yang Exp $
    % No spike found in vm6c_b1.4dfp.img
    % rec vm6c_b1.4dfp.img  Thu May 18 17:15:18 2000  avi
    % /data/petsun4/data1/solaris/imato4dfp2 -fy /data/petsun23/vm6c/siem_im/bold1/5250 7 7 vm6c_b1 
    % $Id: imato4dfp2.c,v 1.12 2000/05/05 00:56:18 avi Exp $
    % patient_id:		vm6c
    % institution:		Washington University
    % manufacturer_model:	MAGNETOM VISION
    % parameter_file_name:	Initialized by sequence
    % sequence_file_name:	/usr/users/tec/nbea_uc_tg2.ekc
    % sequence_description:	ep_fid   90	TR    135.2	TE   37.0/1
    % tilts:			Cor>Tra -12	           
    % 4dfp_dimensions:	64        64        18        128       
    % voxel_dimensions:	3.000000  3.000000  3.000000  
    % scan_date:		22-FEB-1999
    % scan_time:		14:06:33-14:06:33
    % endrec Thu May 18 17:15:18 2000  avi
    % endrec
    % endrec Thu May 18 17:16:26 2000  avi
    % ----------------------------------------------------------------------------------------------------
    % 
    % The following is vm6c_b1_rmsp_dbnd.4dfp.img.rec filtered through brec.
    % ----------------------------------------------------------------------------------------------------
    % 1rec vm6c_b1_rmsp_dbnd.4dfp.img  Thu May 18 17:16:23 2000  avi
    % 1      /data/petsun4/data1/solaris/deband_4dfp -n4 vm6c_b1_rmsp 
    % 1      $Id: deband_4dfp.c,v 1.8 1999/11/20 00:55:49 avi Exp $
    % 1      Frame          1 slice multipliers: even=0.837060 odd=1.162940
    % 1      Frame          2 slice multipliers: even=0.997099 odd=1.002901
    % 1      Frame          3 slice multipliers: even=0.985484 odd=1.014516
    % 1      Frame          4 slice multipliers: even=0.986583 odd=1.013417
    % 1      Functional frame slice multipliers: even=0.986982 odd=1.013018
    % 2      rec vm6c_b1_rmsp.4dfp.img  Thu May 18 17:16:13 2000 avi
    % 2            /data/petsun4/data1/solaris/rmspike_4dfp -n4 -x33 vm6c_b1 
    % 2            $Header: /data/petsun4/src_solaris/rmspike_4dfp/RCS/rmspike_4dfp.c,v 2.6 1997/05/23 00:49:24 yan
    % 2            No spike found in vm6c_b1.4dfp.img
    % 3            rec vm6c_b1.4dfp.img  Thu May 18 17:15:18 2000  avi
    % 3                  /data/petsun4/data1/solaris/imato4dfp2 -fy /data/petsun23/vm6c/siem_im/bold1/5250 7 7 vm6c
    % 3                  $Id: imato4dfp2.c,v 1.12 2000/05/05 00:56:18 avi Exp $
    % 3                  patient_id:             vm6c
    % 3                  institution:            Washington University
    % 3                  manufacturer_model:     MAGNETOM VISION
    % 3                  parameter_file_name:    Initialized by sequence
    % 3                  sequence_file_name:     /usr/users/tec/nbea_uc_tg2.ekc
    % 3                  sequence_description:   ep_fid   90     TR    135.2     TE   37.0/1
    % 3                  tilts:                  Cor>Tra -12                
    % 3                  4dfp_dimensions:        64        64        18        128       
    % 3                  voxel_dimensions:       3.000000  3.000000  3.000000  
    % 3                  scan_date:              22-FEB-1999
    % 3                  scan_time:              14:06:33-14:06:33
    % 3            endrec Thu May 18 17:15:18 2000  avi
    % 2      endrec
    % 1endrec Thu May 18 17:16:26 2000  avi
    % ----------------------------------------------------------------------------------------------------

	%  $Revision$
 	%  was created 30-Apr-2018 16:37:44 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
    properties (Constant) 
        FILETYPE      = '4DFP'
        FILETYPE_EXT  = '.4dfp.hdr'
        FOURDFP_EXT   = '.4dfp.hdr';
        SUPPORTED_EXT = {'.4dfp.hdr' '.4dfp.ifh' '.4dfp.img'};
    end
    
    properties (Dependent)
        ifh
        ifh_orientation
        imgrec
    end
    
    methods (Static)
        function e       = defaultFilesuffix
            e =  mlfourdfp.FourdfpInfo.FILETYPE_EXT;
        end
        function [X,hdr] = exportFourdfp(X, hdr)
            %% EXPORTFOURDFP
            %  Use to maintain interoperability with output of niftigz_4dfp -4 <in.nii.gz> <out.4dfp.hdr> -N
            %  niftigz_4dfp is not compliant with NIfTI qfac.
            
            import mlfourdfp.FourdfpInfo.*;
            if (hdrIsReasonableSurfer(hdr))
                switch ndims(X)
                    case 3
                        X = permute(X, [1 3 2]); % rl, pa, si with respect to fsleyes voxel/world orientations
                    case 4
                        X = permute(X, [1 3 2 4]); % rl, pa, si with respect to fsleyes voxel/world orientations
                    otherwise
                end
                %X = flip(X,1);
                X = flip(X,3);
                
                %X = flip(X,1);
                X = flip(X,2);
            else
                X = flip(X,1);
                X = flip(X,2);
            end
            hdr = adjustHdrForExport(hdr);
        end
        function [X,hdr] = exportFourdfpToFreesurferSpace(X, hdr)
            
            import mlfourdfp.FourdfpInfo.*;
            if (hdrIsReasonableSurfer(hdr))
                X = flip(X,3);
                X = flip(X,1);
                X = permute(X, [1 3 2]);
            end
        end
        function X       = exportFourdfpToNIfTI(X, varargin)
        end
        function [X,hdr] = exportFreeSurferSpaceToFourdfp(X, hdr)
            %% EXPORTFREESURFERSPACETOFOURDFP
            %  Use to maintain interoperability with output of niftigz_4dfp -4 <in.nii.gz> <out.4dfp.hdr> -N
            %  niftigz_4dfp is not compliant with NIfTI qfac.
                  
            import mlfourdfp.FourdfpInfo.*;
            if (hdrIsReasonableSurfer(hdr))
                X = permute(X, [1 3 2]); % rl, pa, si with respect to fsleyes voxel/world orientations
                X = flip(X,1);
                X = flip(X,3);
            end
            hdr = adjustHdrForExport(hdr);    
            
            % eigen flips
            % X = flip(X,1); % rl, pa, si -> LR, pa, si
            % X = flip(X,2); % rl, pa, si -> rl, AP, si 
            % X = flip(X,3); % rl, pa, si -> rl, pa, IS
            % rl, AP, IS is target for 4dfp
        end
        function tf      = hdrIsReasonableNii(~)
            tf = false;
        end
        function tf      = hdrIsReasonableSurfer(hdr)
            tf = all(abs(hdr.hist.srow_x( 2:3 ))             < 0.05) && ...
                 all(abs(hdr.hist.srow_y( 1:2 ))             < 0.05) && ...
                 all(abs(hdr.hist.srow_z([true false true])) < 0.05);
        end
        function X       = importFourdfp(X, varargin)
            ip = inputParser;
            addOptional(ip, 'orientation', 2, @isnumeric);
            parse(ip, varargin{:});
                      
            switch (ip.Results.orientation)
                case 2 % transverse
                    X = flip(X,2);
                case 3 % coronal
                    X = flip(flip(X,2),3);
                case 4 % sagittal
                    X = flip(flip(flip(X,1),2),3);
                otherwise
                    error('mlfourd:unsupportedSwitchcase', 'NIfTId.flip_nii');
            end
        end
    end

	methods 
        
        %% GET, SET
        
        function g    = get.ifh(this)
            g = this.ifh_;
        end 
        function g    = get.ifh_orientation(this)
            if (lstrfind(lower(this.Orientation), 'transverse'))
                g = 2;
                return
            end
            if (lstrfind(lower(this.Orientation), 'coronal'))
                g = 3;
                return
            end
            if (lstrfind(lower(this.Orientation), 'sagittal'))
                g = 4;
                return
            end
            
            error('mlfourdfp:ValueError', 'FourdfpInfo.get.ifh_orientation');
        end 
        function this = set.ifh(this, s)
            this.ifh_ = s;
        end 
        function g    = get.imgrec(this)
            g = this.imgrec_;
        end 
        function this = set.imgrec(this, s)
            this.imgrec_ = s;
        end 
        
        %%
                
        function fqfn = fqfileprefix_4dfp_hdr(this)
            fqfn = [this.fqfileprefix '.4dfp.hdr'];
        end
        function fqfn = fqfileprefix_4dfp_ifh(this)
            fqfn = [this.fqfileprefix '.4dfp.ifh'];
        end
        function fqfn = fqfileprefix_4dfp_img(this)
            fqfn = [this.fqfileprefix '.4dfp.img'];
        end
        function fqfn = fqfileprefix_4dfp_imgrec(this)
            fqfn = [this.fqfileprefix '.4dfp.img.rec'];
        end
        
        function [X,untouch,hdr] = analyze75read(this)
            %% calls mlniftitools.load_nii
            
            jimmy = this.load_nii; % struct
            X = jimmy.img; % identical to nii0 = load_untouch_nii(<from nifti_4dfp -n>); nii0.img
            X = flip(X, 1); % storage order := Neurological, removing Analyze conventions  
            X = this.importFourdfp(X, struct(this.ifh).orientation);
            X = this.ensureDatatype(X, this.datatype_);
            untouch = false;
            hdr = this.adjustHdr(jimmy.hdr);
        end
        function hdr  = recalculateHdrHistOriginator(this, hdr)
            [~,srow] = this.imgrecread;
            if (~isempty(srow) && ~this.N)
                hdr.hist.originator(1) = -srow(1,4);
                hdr.hist.originator(2) = -srow(2,4);
                hdr.hist.originator(3) = -srow(3,4);
            end
            
            hdr.dime.xyzt_units = 2+8; % mm, sec; see also mlniftitools.extra_nii_hdr
            hdr.hist.qform_code = 0;
            hdr.hist.sform_code = 1;
            
                                    % a = 0.5  * sqrt(1 + trace(R)) = 0.5 + sqrt(2);
            hdr.hist.quatern_b = 0; % 0.25 * (R(3,2) - R(2,3)) / a;
            hdr.hist.quatern_c = 0; % 0.25 * (R(1,3) - R(3,1)) / a;
            hdr.hist.quatern_d = 0; % 0.25 * (R(2,1) - R(1,2)) / a;
            
            if (isfield(hdr.hist, 'originator'))
                hdr.hist.qoffset_x = (1-hdr.hist.originator(1))*hdr.dime.pixdim(2);
                hdr.hist.qoffset_y = (1-hdr.hist.originator(2))*hdr.dime.pixdim(3);
                hdr.hist.qoffset_z = (1-hdr.hist.originator(3))*hdr.dime.pixdim(4);            

                % for compliance with NIfTI format
                srow = [[hdr.dime.pixdim(2) 0 0           ( 1-hdr.hist.originator(1))*hdr.dime.pixdim(2)]; ...
                        [0 hdr.dime.pixdim(3) 0           ( 1-hdr.hist.originator(2))*hdr.dime.pixdim(3)]; ...
                        [0 0 hdr.dime.pixdim(4)*this.qfac ( 1-hdr.hist.originator(3))*hdr.dime.pixdim(4)]];
                hdr.hist.srow_x = srow(1,:);
                hdr.hist.srow_y = srow(2,:);
                hdr.hist.srow_z = srow(3,:); 
            end
                       
            hdr = mlniftitools.extra_nii_hdr(hdr);
            hdr = this.adjustDime(hdr);
            hdr = this.adjustHist(hdr);
        end        
        function save(this)
            assert(lstrfind(this.filesuffix, '.4dfp'));
            try                
                warning('off', 'MATLAB:structOnObject');            
                this.ifh.fqfileprefix = this.fqfileprefix;
                this.ifh.save(mlfourd.ImagingFormatContext(this));
                this.imgrec.fqfileprefix = this.fqfileprefix;
                this.imgrec.add('mlfourdfp.FourdfpInfo.save %s', this.fqfileprefix);
                this.imgrec.save;
                warning('on', 'MATLAB:structOnObject');
            catch ME
                dispexcept(ME, ...
                    'mlfourdfp:IOError', ...
                    'FourdfpInfo.save erred while attempting to save %s', this.fqfilename);
            end
        end
        
 		function this = FourdfpInfo(varargin)
 			%% FOURDFPINFO
 			%  @param filename is required.
 			
            this = this@mlfourd.Analyze75Info(varargin{:});
            
            this.raw_.dim = this.hdr_.dime.dim;            
            this.raw_.pixdim = this.hdr_.dime.pixdim;
            this.ifh_ = this.initialIfh;
            this.imgrec_ = this.initialImgRec;
 		end
    end 
    
    %% PROTECTED
    
    methods (Access = protected)
    end
    
    %% PRIVATE
    
    properties (Access = private)
        ifh_ % struct
        imgrec_ % cell
    end
    
    methods (Static, Access = private)
        function hdr = adjustHdrForExport(hdr)
            %% ADJUSTHDRFOREXPORT
            %  Use to maintain interoperability with output of niftigz_4dfp -4 <in.nii.gz> <out.4dfp.hdr> -N
            %  niftigz_4dfp is not compliant with NIfTI qfac.
            
            if (mlfourdfp.FourdfpInfo.hdrIsReasonableNii(hdr))
                return
            end
            
            hdr.hist.qform_code = 0;
            hdr.hist.sform_code = 1;
            
            fx = 1;
            fy = 1;
            fz = 1;
            %fx = hdr.dime.pixdim(2);
            %fy = hdr.dime.pixdim(3);
            %fz = hdr.dime.pixdim(4);
            
            assert(mlpipeline.ResourcesRegistry.instance().defaultN);
            
            srow = [ [hdr.dime.pixdim(2) 0 0 (1-hdr.hist.originator(1))*fx]; ...
                     [0 hdr.dime.pixdim(3) 0 (1-hdr.hist.originator(2))*fy]; ...
                     [0 0 hdr.dime.pixdim(4) (1-hdr.hist.originator(3))*fz] ];
            
            hdr.hist.srow_x = -srow(1,:);
            hdr.hist.srow_y =  srow(2,:);
            hdr.hist.srow_z =  srow(3,:);
        end
    end
    
    methods (Access = private)
        function tf   = consistentByteOrder(this, bo)
            if (strcmp(this.ByteOrder, 'ieee-le'))
                tf = strcmp(bo, 'littleendian');
            end
            if (strcmp(this.ByteOrder, 'ieee-be'))
                tf = strcmp(bo, 'bigendian');
            end
        end
        function tf   = consistentOrientation(this, o)
            switch (lower(this.Orientation(1:7)))
                case 'transve'
                    tf = (2 == o);
                case 'coronal'
                    tf = (3 == o);
                case 'sagitta'
                    tf = (4 == o);
                otherwise
                    error('mlfourdfp:unsupportedSwitchcase', 'FourdfpInfo.assertOrientation');
            end
        end    
        function p    = ifhread(this)
            import mlfourdfp.*;
            p = IfhParser.load([this.fqfileprefix '.4dfp.ifh'], 'N', this.N);
            s = struct(p);
            
            %% possibly breaking pipeline
%             assert(lstrfind(lower(this.ImgDataType), lower(s.number_format)), ...
%                 'mlfourdfp:inconsistentParamValue', ...
%                 'FourdfpInfo.ifhread.ifh.number_format->%s != %s', ...
%                 s.number_format, this.ImgDataType);
            
%            [~,fp] = myfileparts(s.name_of_data_file);
%            assert(lstrfind(fp, this.fileprefix), ...
%               'mlfourdfp:inconsistentParamValue', ...
%               'FourdfpInfo.ifhread.ifh.name_of_data_file->%s ~= %s', ... 
%                s.name_of_data_file, this.fileprefix);
%            assert(s.number_of_bytes_per_pixel == this.BitDepth/8, ...
%                'mlfourdfp:inconsistentParamValue', ...
%                'FourdfpInfo.ifhread.ifh.number_of_bytes_per_pixel->%g != %g', ...
%                 s.number_of_bytes_per_pixel, this.BitDepth/8);
%            assert(this.consistentByteOrder(s.imagedata_byte_order), ...
%                'mlfourdfp:inconsistentParamValue', ...
%                'FourdfpInfo.ifhread.ifh.imagedata_byte_order->%s', s.imagedata_byte_order);
            assert(this.consistentOrientation(s.orientation), ...
                'mlfourdfp:inconsistentParamValue', ...
                'FourdfpInfo.ifhread.ifh.orientation->%g', s.orientation);   
        end
        function [cc,sr] = imgrecread(this)
            p = mlfourdfp.ImgRecParser.load([this.fqfileprefix '.4dfp.img.rec']);
            cc = p.cellContents;
            sr = p.commonSform;
        end
        function ii   = initialIfh(this)
            if (~lexist([this.fqfileprefix '.4dfp.ifh'], 'file'))
                ii = mlfourdfp.IfhParser.constructDenovo( ...
                    this.hdr, ...
                    'fqfileprefix', this.fqfileprefix, ...
                    'orientation', 2, ...
                    'N', this.N);
                return
            end            
            ii = this.ifhread;
        end
        function irl  = initialImgRec(this)
            irl = mlfourdfp.ImgRecLogger(this.fqfileprefix, this);
        end
    end
    
    %% 
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

