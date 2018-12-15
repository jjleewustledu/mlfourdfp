classdef InnerFourdfp < handle & mlfourd.AbstractInnerImagingFormat
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
        
        function x = get.hdxml(~)
            x = '';
        end
        function g = get.ifh(this)
            if (~isprop(this.imagingInfo_, 'ifh'))
                g = [];
                return
            end
            g = this.imagingInfo_.ifh;
        end
        function     set.ifh(this, s)
            if (~isprop(this.imagingInfo_, 'ifh'))
                return
            end
            assert(isa(s, mlfourdfp.IfhParser));
            this.imagingInfo_.ifh = s;
        end
        function o = get.orient(this)
            if (~isempty(this.orient_))
                o = this.orient_;
                return
            end
            if (lexist(this.fqfilename, 'file'))
                [~, o] = mlbash(['fslorient -getorient ' this.fqfileprefix '.4dfp.hdr']);
                o = strtrim(o);
                return
            end
            o = '';
        end
        function u = get.untouch(this)
            u = this.imagingInfo_.untouch;
        end
        function     set.untouch(this, s)
            this.imagingInfo_.untouch = logical(s);
        end
        
        %%
        
        function deleteExisting(~, fn)
            deleteExisting_4dfp(fn);
        end
        function this = mutateInnerImagingFormatByFilesuffix(this)
            import mlfourd.* mlfourdfp.* mlsurfer.*;  
            hdr_ = this.hdr;
            switch (this.filesuffix)
                case FourdfpInfo.SUPPORTED_EXT
                case [NIfTIInfo.SUPPORTED_EXT]
                    deleteExisting([this.fqfileprefix '.nii*']);
                    this.img = FourdfpInfo.exportFourdfpToNIfTI(this.img, struct(this.ifh).orientation);
                    info = NIfTIInfo(this.fqfilename, ...
                        'datatype', this.datatype, 'ext', this.imagingInfo.ext, 'filetype', this.imagingInfo.filetype, 'N', this.N , 'untouch', false, 'hdr', this.hdr);
                    this = InnerNIfTI(info, ...
                       'creationDate', this.creationDate, 'img', this.img, 'label', this.label, 'logger', this.logger, ...
                       'orient', this.orient_, 'originalType', this.originalType, 'seriesNumber', this.seriesNumber, ...
                       'separator', this.separator, 'stack', this.stack, 'viewer', this.viewer);
                    this.imagingInfo.hdr = hdr_;
                case MGHInfo.SUPPORTED_EXT 
                    deleteExisting([this.fqfileprefix '.mgz']);
                    deleteExisting([this.fqfileprefix '.mgh']);
                    [this.img_,hdr_] = FourdfpInfo.exportFourdfpToFreesurferSpace(this.img_, hdr_);
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
        
        function this = InnerFourdfp(varargin)
 			%  @param imagingInfo is an mlfourd.ImagingInfo object and is required; it may be an aufbau object.
            
            this = this@mlfourd.AbstractInnerImagingFormat(varargin{:});
            this = this.adjustQOffsets;
            this = this.adjustSRows;
        end
    end 
    
    %% PROTECTED
    
    methods (Access = protected)
        function [s,r] = viewExternally(this, app, varargin)
            s = []; r = '';
            that = copy(this); % avoid side effects
            assert(0 == mlbash(sprintf('which %s', app)), ...
                'mlfourdfp:externalAppNotFound', ...
                'InnerFourdfp.viewExternally could not find %s', app);
            tmp = that.tempFqfilename;
            that.fqfilename = tmp;
            try
                that.save; % always save temp; internal img likely has changed from img on filesystem
                v = mlfourdfp.Viewer(app);
                tmp = [tmp varargin];
                [s,r] = v.aview(tmp{:});
                that.deleteExisting(tmp);
            catch ME
                dispexcept(ME, 'mlfourdfp:viewerError', ...
                    'InnerFourdfp.viewExternally called mlbash with %s; \nit returned s->%i, r->%s', ...
                    app, s, r);
            end
        end        
    end

    %% HIDDEN
    
    methods (Hidden) 
        function save__(this)
            that = copy(this);
            assert(lstrfind(this.filesuffix, '.4dfp'));
            [that.img_,hdr_] = mlfourdfp.FourdfpInfo.exportFourdfp(that.img_, that.imagingInfo.hdr); %% KLUDGE
            that.imagingInfo_.hdr = hdr_;
            try                
                warning('off', 'MATLAB:structOnObject');
                ana = mlniftitools.make_ana(single(that.img_), that.mmppix);
                mlniftitools.save_untouch_nii(ana, that.fqfileprefix_4dfp_hdr);                
                that.imagingInfo_.ifh.fqfileprefix = that.fqfileprefix;
                that.imagingInfo_.ifh.save(that);
                that.imagingInfo_.imgrec.fqfileprefix = that.fqfileprefix;
                that.imagingInfo_.imgrec.save;
                warning('on', 'MATLAB:structOnObject');
            catch ME
                dispexcept(ME, ...
                    'mlfourdfp:IOError', ...
                    'InnerFourdfp.saveByNifti4dfp erred while attempting to save %s', that.fqfilename);
            end
        end
    end
    
    %% PRIVATE
    
    methods (Access = private)
        function this = adjustQOffsets(this)
            this.hdr.hist.qoffset_x = this.hdr.hist.qoffset_x / this.hdr.dime.pixdim(2);
            this.hdr.hist.qoffset_y = this.hdr.hist.qoffset_y / this.hdr.dime.pixdim(3);
            this.hdr.hist.qoffset_z = this.hdr.hist.qoffset_z / this.hdr.dime.pixdim(4);
        end
        function this = adjustSRows(this)
            this.hdr.hist.srow_x(4) = this.hdr.hist.srow_x(4) / this.hdr.dime.pixdim(2);
            this.hdr.hist.srow_y(4) = this.hdr.hist.srow_y(4) / this.hdr.dime.pixdim(2);
            this.hdr.hist.srow_z(4) = this.hdr.hist.srow_z(4) / this.hdr.dime.pixdim(2);
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

