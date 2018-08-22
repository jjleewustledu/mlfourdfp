classdef InnerFourdfp < mlfourd.AbstractInnerImagingFormat
	%% INNERFOURDFP  

	%  $Revision$
 	%  was created 21-Jul-2018 23:13:59 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties (Dependent)
 		hdxml
        ifh
        orient % RADIOLOGICAL, NEUROLOGICAL
        untouch
 	end

	methods (Static)
        function this = create(fn, varargin)
            import mlfourdfp.*;
            this = InnerFourdfp( ...
                InnerFourdfp.createImagingInfo(fn, varargin{:}), varargin{:});
        end
        function info = createImagingInfo(fn, varargin)
            info = mlfourdfp.FourdfpInfo(fn, varargin{:});
        end
        function [s,finfo] = imagingInfo2struct(fn, varargin)
            fn = [myfileprefix(fn) '.4dfp.hdr'];
            finfo = mlfourdfp.FourdfpInfo(fn, varargin{:});
            nii = finfo.make_nii;
            s = struct( ...
                'hdr', nii.hdr, ...
                'filetype', 2, ...
                'fileprefix', fn, ...
                'machine', finfo.machine, ...
                'ext', [], ...
                'img', nii.img, ...
                'untouch', nii.untouch);  
            
            %% DEBUGGING
            %save_nii(nii, 'test.nii.gz');
            %system('fsleyes test.nii.gz');
        end
    end
    
    methods
        
        %% GET/SET
        
        function x    = get.hdxml(~)
            x = '';
        end
        function g    = get.ifh(this)
            if (~isprop(this.imagingInfo_, 'ifh'))
                g = [];
                return
            end
            g = this.imagingInfo_.ifh;
        end
        function this = set.ifh(this, s)
            if (~isprop(this.imagingInfo_, 'ifh'))
                return
            end
            assert(isa(s, mlfourdfp.IfhParser));
            this.imagingInfo_.ifh = s;
        end
        function o    = get.orient(this)
            if (~isempty(this.orient_))
                o = this.orient_;
                return
            end
            if (lexist(this.fqfilename, 'file') && lstrfind(this.filesuffix, '.hdr'))
                [~, o] = mlbash(['fslorient -getorient ' this.fqfileprefix]);
                o = strtrim(o);
                return
            end
            o = '';
        end
        function u    = get.untouch(this)
            u = this.imagingInfo_.untouch;
        end
        function this = set.untouch(this, s)
            this.imagingInfo_.untouch = logical(s);
        end
        
        %%
        
        function deleteExisting(~, fn)
            deleteExisting_4dfp(fn);
        end
        
        function this = InnerFourdfp(varargin)
 			%  @param imagingInfo is an mlfourd.ImagingInfo object and is required; it may be an aufbau object.
            
            this = this@mlfourd.AbstractInnerImagingFormat(varargin{:});
        end
    end 
    
    %% PROTECTED
    
    methods (Access = protected)
        function this = mutateInnerImagingFormatByFilesuffix(this)
            import mlfourd.* mlfourdfp.* mlsurfer.*;  
            hdr_ = this.hdr;
            switch (this.filesuffix)
                case FourdfpInfo.SUPPORTED_EXT
                    [this.img,hdr_] = FourdfpInfo.exportFourdfp(this.img, hdr_);
                    this.imagingInfo.hdr = hdr_;
                case [NIfTIInfo.SUPPORTED_EXT '.hdr']
                    this.img = FourdfpInfo.exportFourdfpToNIfTI(this.img, this.ifh.asstruct.orientation);
                    info = NIfTIInfo(this.fqfilename, ...
                        'datatype', this.datatype, 'ext', this.imagingInfo.ext, 'filetype', this.imagingInfo.filetype, 'N', this.N , 'untouch', false, 'hdr', this.hdr);
                    this = InnerNIfTI(info, ...
                       'creationDate', this.creationDate, 'img', this.img, 'label', this.label, 'logger', this.logger, ...
                       'orient', this.orient_, 'originalType', this.originalType, 'seriesNumber', this.seriesNumber, ...
                       'separator', this.separator, 'stack', this.stack, 'viewer', this.viewer);
                    this.imagingInfo.hdr = hdr_;
                case MGHInfo.SUPPORTED_EXT 
                    this.img = FourdfpInfo.exportFourdfpToNIfTI(this.img, this.ifh.asstruct.orientation);
                    info = MGHInfo(this.fqfilename, ...
                        'datatype', this.datatype, 'ext', this.imagingInfo.ext, 'filetype', this.imagingInfo.filetype, 'N', this.N , 'untouch', false, 'hdr', this.hdr);
                    this = InnerMGH(info, ...
                       'creationDate', this.creationDate, 'img', this.img, 'label', this.label, 'logger', this.logger, ...
                       'orient', this.orient_, 'originalType', this.originalType, 'seriesNumber', this.seriesNumber, ...
                       'separator', this.separator, 'stack', this.stack, 'viewer', this.viewer);
                    this.imagingInfo.hdr = hdr_;
                otherwise
                    error('mlfourd:unsupportedSwitchcase', ...
                        'InnerNIfTI.filesuffix->%s', this.filesuffix);
            end
        end
        function [s,r] = viewExternally(this, app, varargin)
            s = []; r = '';
            try
                assert(0 == mlbash(sprintf('which %s', app)), ...
                    'mlfourdfp:externalAppNotFound', ...
                    'InnerFourdfp.viewExternally could not find %s', app);
                tmp = this.tempFqfilename;
                this.saveas(tmp); % always save temp; internal img likely has changed from img on filesystem
                v = mlfourdfp.Viewer(app);
                [s,r] = v.aview([tmp varargin{:}]);
                this.deleteExisting(tmp);
            catch ME
                handexcept(ME, 'mlfourdfp:viewerError', ...
                    'InnerFourdfp.viewExternally called mlbash with %s; \nit returned s->%i, r->%s', ...
                    app, s, r);
            end
        end        
    end

    %% HIDDEN
    
    methods (Hidden) 
        function save__(this)
            warning('off', 'MATLAB:structOnObject');
            try
                this.img_ = single(this.img_);
                % this.img = flip(this.img, 2);
                mlniftitools.save_nii(struct(this), this.fqfileprefix_4dfp_hdr);
                this.imagingInfo_.ifh.fqfileprefix = this.fqfileprefix;
                this.imagingInfo_.ifh.save(this);
                this.imagingInfo_.imgrec.fqfileprefix = this.fqfileprefix;
                this.imagingInfo_.imgrec.save;
            catch ME
                dispexcept(ME, ...
                    'mlfourdfp:IOError', ...
                    'InnerFourdfp.saveByNifti4dfp erred while attempting to save %s', this.fqfilename);
            end
            warning('on', 'MATLAB:structOnObject');
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

