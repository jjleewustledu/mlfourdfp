classdef FourdfpVisitor 
	%% FOURDFPVISITOR  

	%  $Revision$
 	%  was created 01-Mar-2016 18:43:03
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.307022 (R2016a) Prerelease for MACI64.
 	

	properties (Constant)
 		FOURDFP_HOSTS = {'touch3'}
    end
    
    methods (Static)
        function tf    = lexist(fdfp)
            tf = lexist(fdfp);
            if (~tf)
                tf = lexist([fdfp '.4dfp.hdr']) && ...
                     lexist([fdfp '.4dfp.ifh']) && ...
                     lexist([fdfp '.4dfp.img']);
            end
        end
        function [s,r] = lns(fqfn)
            if (~lexist(basename(fqfn), 'file'))
                [s,r] = mlbash(sprintf('ln -s %s', fqfn));
            end
        end
        function [s,r] = lns_4dfp(fqfp)
            ext = { '.hdr' '.ifh' '.img' '.img.rec' };
            for e = 1:length(ext) 
                if (~lexist([basename(fqfp) '.4dfp' ext{e}], 'file'))
                    [s,r] = mlbash(sprintf('ln  -s %s.4dfp%s', fqfp, ext{e}));
                end
            end
        end
        function [s,r] = mkdir(pth)
            if (~isdir(pth))
                [s,r] = mlbash(sprintf('mkdir -p %s', pth));
            end
        end
        function [s,r] = pushd(pth)
            [s,r] = mlbash(sprintf('pushd %s', pth));
        end
        function [s,r] = popd
            [s,r] = mlbash('popd');
        end
    end

	methods
 		function this = FourdfpVisitor(varargin)
 			%% FOURDFPVISITOR
 			%  Usage:  this = FourdfpVisitor()

 			this.assertPlatform;
        end
        
        function resolveTracer(this, varargin)
            t4_resolver_subject;
            t4_resolver_crossmodal;
            t4_resolver_top;
            t4_resolver;
            crop_4dfp;
            copy_4dfp;
            msktgen_4dfp_mpr;
            msktgen_4dfp_init;
            extract_frame_4dfp;
            frame_reg;
            t4_img_4dfp;
            t4_resolve_and_paste;
            msktgen_4dfp_resolved;
            t4_resolver_teardown;
        end
        
        function [s,r] = actmapf_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'format',      @ischar);
            addRequired( ip, 'input',       @this.lexist);
            addParameter(ip, 'options', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.actmapf_4dfp__(sprintf( ...
                '%s %s %s', ip.Results.format, ip.Results.input, ip.Results.options));
        end
        function [s,r] = extract_frame_4dfp(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fdfp',        @(x) this.lexist(x));
            addRequired(ip, 'frame',       @isnumeric);
            addOptional(ip, 'options', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.extract_frame_4dfp__( ...
                sprintf('%s %i %s', ip.Results.fdfp, ip.Results.frame, ip.Results.options));
        end
        function [s,r] = gauss_4dfp(this, varargin)
            ip = inputParser;
            addRequired(ip, 'input',            @this.lexist);
            addRequired(ip, 'f_half',           @isnumeric);
            outroot = ...
                sprintf('%s_g%s', varargin{1}, num2str(varargin{2}));
            addOptional(ip, 'outroot', outroot, @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.gauss_4dfp__(sprintf('%s %g %s', ip.Results.input, ip.Results.f_half, ip.Results.outroot));
        end
        function [s,r] = imgblur_4dfp(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fdfp',        @this.lexist);
            addRequired(ip, 'fwhm',        @isnumeric);
            addOptional(ip, 'options', '', @ischar);            
            parse(ip, varargin{:});
            
            [s,r] = this.imgblur_4dfp__( ...
                sprintf('%s %s %g', ip.Results.options, ip.Results.fdfp, ip.Results.fwhm));
        end
        function [s,r] = imgreg_4dfp(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fdfp0', @this.lexist);
            addRequired(ip, 'mask0', @ischar);
            addRequired(ip, 'fdfp1', @this.lexist);
            addRequired(ip, 'mask1', @ischar);
            addRequired(ip, 't4',    @ischar);
            addRequired(ip, 'mode',  @isnumeric);
            addOptional(ip, 'log', 'null', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.imgreg_4dfp__(sprintf('%s %s %s %s %s %i >> %s', ...
                ip.Results.fdfp0, ...
                ip.Results.mask0, ...
                ip.Results.fdfp1, ...
                ip.Results.mask1, ...
                ip.Results.t4, ...
                ip.Results.mode, ...
                ip.Results.log));
        end
        function [s,r] = mpr2atl_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'in',             @(x) lexist(x, 'file') || lexist([x '.4dfp.img'], 'file'));
            addParameter(ip, 'options',   '',  @ischar);
            addParameter(ip, 'log',       '',  @ischar);
            parse(ip, varargin{:});
            
            cmd = '%s %s';
            if (~isempty(ip.Results.log))
                cmd = [cmd ' &> ' ip.Results.log]; 
            end
            [s,r] = this.mpr2atl_4dfp__(sprintf( ...
                cmd, ip.Results.in, ip.Results.options));
        end
        function [s,r] = msktgen_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'in',             @(x) lexist(x, 'file') || lexist([x '.4dfp.img'], 'file'));
            addOptional( ip, 'threshold', 200, @isnumeric);
            addParameter(ip, 'options',   '',  @ischar);
            addParameter(ip, 'log',       '',  @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.msktgen_4dfp__( ...
                sprintf('%s %i %s %s', ip.Results.in, ip.Results.threshold, ip.Results.options, this.log(ip.Results)));
        end
        function [s,r] = paste_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'inlist',      @lexist);
            addRequired( ip, 'outfile',     @ischar);
            addParameter(ip, 'options', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.paste_4dfp__( ...
                sprintf(' %s %s %s', ip.Results.options, ip.Results.inlist, ip.Results.outfile));
        end
        function [s,r] = t4img_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 't4',          @this.lexist);
            addRequired( ip, 'in',          @this.lexist);
            addRequired( ip, 'out',         @ischar);
            addParameter(ip, 'options', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.t4img_4dfp__( ...
                sprintf('%s %s %s %s', ip.Results.t4, ip.Results.in, ip.Results.out, ip.Results.options));
        end
        function [s,r] = t4_inv(this, varargin)
            ip = inputParser;
            addRequired(ip, 'in',  @(x) lexist(x, 'file'));
            addRequired(ip, 'out', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.t4_inv__(sprintf('%s %s', ip.Results.in, ip.Results.out));
            try
                dbbash(sprintf('chmod 664 %s', ip.Results.out));
            catch ME
                handexcept(ME);
            end
        end
        function [s,r] = t4_mul(this, varargin)
            ip = inputParser;
            addRequired(ip, 'A2B', @(x) lexist(x, 'file'));
            addRequired(ip, 'B2C', @(x) lexist(x, 'file'));
            addOptional(ip, 'A2C', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.t4_mul__(sprintf('%s %s %s', ip.Results.A2B, ip.Results.B2C, ip.Results.A2C));
        end
        function [s,r] = t4_resolve(this, varargin)
            ip = inputParser;
            addRequired( ip, 'output',        @ischar);
            addRequired( ip, 'filenames',     @ischar);
            addParameter(ip, 'options',   '', @ischar);
            addParameter(ip, 'log',           @ischar);
            parse(ip, varargin{:});            
            
            [s,r] = this.t4_resolve__( ...
                sprintf(' %s -o%s %s %s', ip.Results.options, ip.Results.output, ip.Results.filenames, this.log(ip.Results)));
        end
    end 
    
    %% PROTECTED
    
    methods (Access = protected)
        function assertPlatform(this)
            if (logical(getenv('DEBUG')))
                return
            end            
            assert(strcmp(computer, 'GLNXA64'));
            [~,hn] = mlbash('hostname');
            assert(lstrfind(hn, this.FOURDFP_HOSTS));
        end        
        function fn = log(~, results)
            assert(isfield(results, 'log'));
            if (~isempty(results.log))
                fn = [' &> ' results.log]; 
                return
            end
            fn = '';            
        end
        
        function [s,r] = actmapf_4dfp__(~, args)
            %% ACTMAPF_4DFP__
            % $Id: actmapf_4dfp.c,v 1.34 2016/02/12 00:59:43 avi Exp $
            % Usage:	actmapf_4dfp <format|fmtfile> <4dfp|conc input>
            %  e.g.,	actmapf_4dfp -zu "3x3(11+4x15-)" b1_rmsp_dbnd_xr3d_norm
            %  e.g.,	actmapf_4dfp -aanatomy -c10 -u "+" ball_dbnd_xr3d.conc
            %  e.g.,	actmapf_4dfp -zu "4x124+" b1_rmsp_dbnd_xr3d -wweights.txt
            % 	        option
            % 	        -a<str>	specify 4dfp output root trailer (default = "actmap")
            % 	        -c<flt>	scale output by specified factor
            % 	        -u	scale weights to unit variance
            % 	        -z	adjust weights to zero sum
            % 	        -R	compute relative modulation (default absolute)
            % 	        -w<weight file>	read (text) weights from specified filename
            % 	        -@<b|l>	output big or little endian (default input endian)
            % N.B.:	    conc files must have extension "conc"
            % N.B.:	    when using weight files 'x' frames in format are not counted
            % N.B.:	    relative modulation images are zeroed where mean intensity < 0.5*whole_image_mode
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('actmapf_4dfp %s', args));
        end
        function [s,r] = crop_4dfp__(~, args)
            %% CROP_4DFP__
            % $Id: crop_4dfp.c,v 1.10 2010/11/12 05:54:27 avi Exp $
            % Usage:	crop_4dfp <(4dfp) inroot> [(4dfp) outroot]
            % 	        option
            % 	        -<x|y|z><int>[to[<int>]	specify x y z crop limits
            % 	        -s<x|y|z><int>	scroll specified axis by specified number of pixels
            % 	        -f	interpret specifications under 4dfp<->analyze flips
            % 	        -Z	zero voxels instead of physically cropping
            % 	        -@<b|l>	output big or little endian (default input endian)
            % N.B.:	    crop limit indices count from 1
            % N.B.:	    scrolling is done after cropping
            % N.B.:	    default (4dfp) output root is <(4dfp) inroot>"_crop"
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('crop_4dfp %s', args));
        end
        function [s,r] = extract_frame_4dfp__(~, args)
            %% EXTRACT_FRAME_4DFP__
            % usage:	extract_frame_4dfp <(4dfp) stack> <(int) frame>
            % e.g.:	    extract_frame_4dfp CDR.5to1+ 3
            % 	        option
            % 	        -o<str>	specifiy output 4dfp fileroot (default = <stack>_frame<frame>)
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('extract_frame_4dfp %s', args));
        end
        function [s,r] = gauss_4dfp__(~, args)
            %% GAUSS_4DFP__
            % $Id: gauss_4dfp.c,v 1.16 2007/05/04 01:26:56 avi Exp $
            % Usage:	gauss_4dfp <4dfp|conc input> f_half [outroot]
            %  e.g.,	gauss_4dfp pt349_study9to9 0.1
            %  e.g.,	gauss_4dfp p1234ho5 0.7 p1234ho5_g7
            % 	        options
            % 	        -@<b|l>	output big or little endian (default input endian)
            % 	        -w	(wrap) suppress x and y padding
            % 	        -d	differentiate
            % N.B.:	    f_half is half frequency in 1/cm
            % N.B.:	    default output root is <inroot>_g<10*f_half>
            % N.B.:	    FWHM*f_half = (2ln2/pi) = 0.4412712
            % N.B.:	    conc files must have extension "conc"
            % N.B.:	    user outroot specification not possible with conc files
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('gauss_4dfp %s', args));
        end
        function [s,r] = imgblur_4dfp__(~, args)
            %% IMGBLUR_4DFP__
            % $Id: imgblur_4dfp.c,v 1.9 2009/07/18 01:30:43 avi Exp $
            % Usage:	imgblur_4dfp [options] <image_file> <FWHM_in_mm>
            %  e.g.,	imgblur_4dfp -yz vc345 5.5
            % 	        option
            % 	        -x	slective x blur
            % 	        -y	slective y blur
            % 	        -z	slective z blur
            % 	        -@<b|l>	output big or little endian (default input endian)
            % N.B.:	    default blur is 3D isotropic
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('imgblur_4dfp %s', args));
        end
        function [s,r] = imgreg_4dfp__(~, args)
            %% IMGREG_4DFP__
            % $Id: imgreg_4dfp.c,v 1.15 2007/08/10 03:27:30 avi Exp $
            % Usage:	imgreg_4dfp target_imag target_mask source_imag source_mask t4file mode
            % or:	    imgreg_4dfp target_imag        none source_imag source_mask t4file mode
            % or:	    imgreg_4dfp target_imag        none source_imag        none t4file mode
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('imgreg_4dfp %s', args));
        end
        function [s,r] = mpr2atl_4dfp__(~, args)
            %% MPR2ATL_4DFP__
            % $Id: mpr2atl_4dfp,v 1.27 2014/02/16 02:46:17 avi Exp $
            % Usage:	mpr2atl_4dfp <mpr_anat> [options]
            % e.g.,	    mpr2atl_4dfp vc1234_654-3[.4dfp.img]
            % e.g.,	    mpr2atl_4dfp vc1234_654-3[.4dfp.img] -T/data/petsun23/data1/atlas/NP345_111[.4dfp.img] -S711-2B
            % 	        options
            % 	        711-2<C|O|Y|K|L|G|H|V|F>	specify 711-2? series atlas representative target image
            % 	        -T<target including path>	specify arbitrary     atlas representative target image
            % 	        -S<atlas space>			specify atlas space (default=711-2B space)
            % 	        crossmodal	use cross-modal mpr->target registration
            % 	        useold		suppress recomputation  of existing t4 file
            % 	        redo		suppress initialization of existing t4 file
            % 	        setecho		set echo
            % N.B.:	    <mpr_anat> may include a path, e.g., /data/petmr1/data7/stem9/scout/654-3
            % N.B.:	    <mpr_anat> must be in either ANALYZE short int or 4dfp format; ANALYZE will be converted to 4dfp
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('mpr2atl_4dfp %s', args));
        end
        function [s,r] = msktgen_4dfp__(~, args)
            %% MSKTGEN_4DFP__
            % $Id: msktgen_4dfp,v 1.21 2010/07/27 01:52:58 avi Exp $
            % Usage:	msktgen_4dfp <(4dfp) image> [threshold] -T<target including path>  -S<atlas space>
            % e.g.,	    msktgen_4dfp 4859-5_mpr
            % e.g.,	    msktgen_4dfp 4859-5_mpr -T/data/petsun29/data1/atlas/NP345_111[.4dfp.img] -S711-2B
            % N.B.:	    default threshold = 200
            % 	        Specifiy a higher threshold, e.g., 400 for a tighter mask and vice versa
            % N.B.:	    msktgen_4dfp uses the first legitimate atlas transform t4 file it sees in 
            % 	        the current working directory, i.e., one of <image>_to_711-2*_t4
            %  	        or  one of <image>_to_<target>_t4
            % N.B.:	-S specifies the atlas space to use. Atlas supported currently is 711-2B. -S must be used with -T option
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('msktgen_4dfp %s', args));
        end
        function [s,r] = paste_4dfp__(~, args)
            %% PASTE_4DFP__
            % $Id: paste_4dfp.c,v 1.9 2008/10/17 02:53:48 avi Exp $
            % Usage:	paste_4dfp <inlist> <outfile>
            % 	        option
            % 	        -a	append successive epochs (default average)
            % 	        -p<int>	specify period in frames (default=1)
            % 	        -@<b|l>	output big or little endian (default initial input endian)
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('paste_4dfp %s', args));
        end
        function [s,r] = t4_inv__(~, args)
            %% T4_INV__
            % $Id: t4_inv.c,v 1.1 2007/05/01 01:20:20 avi Exp $
            % Usage:	t4_inv <t4file> [inv_t4file]
            % e.g.,	    t4_inv vm11b_anat_ave_to_vm11b_234-3_t4 [vm11b_234-3_to_vm11b_anat_ave_t4]
            % 	        option
            % 	        -u	suppress (intensity) scale field in output
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('t4_inv %s', args));
        end
        function [s,r] = t4_mul__(~, args)
            %% T4_MUL__
            % $Id: t4_mul.c,v 1.1 2007/04/30 17:37:57 avi Exp $
            % Usage:	t4_mul <left_t4file> <right_t4file> [product_t4file]
            % e.g.:	    t4_mul vm11b_anat_ave_to_vm11b_234-3_t4 vm11b_234-3_to_711-2B_t4 [vm11b_anat_ave_to_711-2B_t4]
            % 	        option
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('t4_mul %s', args));
        end
        function [s,r] = t4img_4dfp__(~, args)
            %% T4IMG_4DFP__
            % $Id: t4img_4dfp,v 1.14 2015/11/03 04:45:16 avi Exp $
            % Usage:	t4img_4dfp <t4file> <imgfile> [outfile]
            %  e.g.,	t4img_4dfp  vce1_mprS_to_711-2B_t4	vce1_mprS.4dfp.img -O222
            %    or,	t4img_4dfp  vce1_mprS_to_711-2B_t4 	vce1_mprS vce_mprS_711-2B -O222
            %    or,	t4img_4dfp  none			vce1_mprS vce1_mprS_222 -O222
            % N.B.:	    4dfp filename extensions are optional
            % N.B.:	    default output filename = <imgfile>t
            % N.B.:	    t4img_4dfp is a wrapper for t4imgs_4dfp; options are listed by t4imgs_4dfp usage
            % N.B.:	    option -n causes fidl ROI names to be copied to the output ifh
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('t4img_4dfp %s', args));
        end
        function [s,r] = t4imgs_4dfp__(~, args)
            %% T4IMGS_4DFP__
            % $Id: t4imgs_4dfp.c,v 1.37 2015/11/03 03:50:32 avi Exp $
            % Usage:	t4imgs_4dfp [options] <inlist> <outfile>
            % 	option
            % 	-z	normalize by sqrt(n) rather than n (for z images)
            % 	-s	interpolate by 3D cubic spline (default is 3D linear)
            % 	-N	output NaN (default 0.0) for undefined values
            % 	-B	internally convert to_711-2A_t4->to_711-2B_t4
            % 	-n	use nearest neighbor interpolation
            % 	-R	suppress creation of rec file
            % 	-O111	output in 111 space instead of default 333.0 space
            % 	-O222	output in 222 space instead of default 333.0 space
            % 	-O333.n	output in 333.n space (y shifted up by n pixels)
            % 	-Omy_image	duplicate dimensions of my_image.4dfp.ifh
            % 	-@<b|l>	output big or little endian (default CPU endian)
            % N.B.:	t4file intensity scale ingnored with option -n
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('t4imgs_4dfp %s', args));
        end
        function [s,r] = t4_resolve__(~, args)
            %% T4_RESOLVE__
            % $Id: t4_resolve.c,v 1.6 2013/09/12 00:31:12 avi Exp $
            % n=0
            % Usage:	t4_resolve <image1> <image2> ...
            % 	        option
            % 	        -v	verbose mode
            % 	        -m	generate mat file output
            % 	        -s	include intensity scale factor in t4 file output
            % 	        -w	weight inversely in proportion to scale in sub file output (sum counts mode)
            % 	        -o<str>	write resolved output with specified fileroot
            % 	        -r<flt>	set VOI rms radius in mm (default=50)
            % N.B.:	    t4_resolve looks for t4 files <image1>_to_<image2>_t4, <image1>_to_<image3>_t4, ...
            % N.B.:	    t4_resolve automatically strips filename extensions when constructing t4 filenames
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('t4_resolve %s', args));
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

