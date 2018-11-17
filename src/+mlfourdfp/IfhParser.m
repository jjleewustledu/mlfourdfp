classdef IfhParser < mlio.AbstractParser
	%% IfhParser parses numerical values to the right or left of a text field-name.
    %  For more fine-grained parsing features, see TextParser.

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.5.0.197613 (R2015a) 
 	%  $Id$   
    
    properties (Constant)
        IFH_EXT = '.4dfp.ifh';
        MIN_VERSION = 3.3
        SUPPORTED_KEYS = { ...
            'version_of_keys' 'number_format' 'conversion_program' 'name_of_data_file' ...
            'number_of_bytes_per_pixel' 'imagedata_byte_order' 'orientation' ...
            'number_of_dimensions' 'matrix_size' 'scaling_factor' 'mmppix' 'center'}
    end
    
	methods (Static)
        function this = constructDenovo(varargin)
            ip = inputParser;
            addRequired( ip, 'hdr', @isstruct);
            addParameter(ip, 'fileprefix', @ischar);
            addParameter(ip, 'orientation', 2, @isnumeric);
            addParameter(ip, 'N', false, @islogical);
            parse(ip, varargin{:});
            hdr = ip.Results.hdr;
            
            assert(isstruct(hdr),        'mlfourdfp:unsupportedInputTypeclass', 'IhfParser.constructDenovo');
            assert(isfield(hdr, 'hk'),   'mlfourdfp:unsupportedInputTypeclass', 'IhfParser.constructDenovo');
            assert(isfield(hdr, 'dime'), 'mlfourdfp:unsupportedInputTypeclass', 'IhfParser.constructDenovo');
            assert(isfield(hdr, 'hist'), 'mlfourdfp:unsupportedInputTypeclass', 'IhfParser.constructDenovo');            
            
            this = mlfourdfp.IfhParser;
            this.fileprefix = ip.Results.fileprefix;
            this.denovo_ = struct( ...
                'version_of_keys', 3.3, ...
                'number_format', 'float', ...
                'conversion_program', 'mlfourdfp.IfhParser.constructDenovo', ...
                'name_of_data_file', ip.Results.fileprefix, ...
                'patient_ID', '', ...
                'date', datestr(now), ...
                'number_of_bytes_per_pixel', 4, ...
                'imagedata_byte_order', this.imagedataByteOrder, ...
                'orientation', ip.Results.orientation, ...
                'number_of_dimensions', 4, ...
                'matrix_size', hdr.dime.dim(2:5), ...
                'global_minimum', hdr.dime.glmin, ...
                'global_maximum', hdr.dime.glmax, ...
                'scaling_factor', hdr.dime.pixdim(2:4), ...
                'slice_thickness', hdr.dime.pixdim(4));
            if (~ip.Results.N), ...
                this.denovo_.mmppix = hdr.dime.pixdim(2:4);
                this.denovo_.center = hdr.hist.originator;     
            end
        end
        function o = imagedataByteOrder
            [~,~,o] = computer;
            if (strcmpi(o, 'L'))
                o = 'littleendian';
            else
                o = 'bigendian';
            end
        end
        function this = load(varargin)
            ip = inputParser;
            addRequired(ip, 'fn', @(x) lexist(x, 'file'));
            addParameter(ip, 'N', mlpet.Resources.instance.defaultN, @islogical);
            parse(ip, varargin{:});
            fn = ip.Results.fn;
            
            [pth, fp, fext] = myfileparts(fn); 
            if (lstrfind(fext, mlfourdfp.IfhParser.IFH_EXT) || ...
                isempty(fext))
                this = mlfourdfp.IfhParser.loadText(fn); 
                this.filepath_   = pth;
                this.fileprefix_ = fp;
                this.filesuffix_ = fext;
                this.N_ = ip.Results.N;
                return 
            end
            error('mlio:unsupportedParam', ...
                'IfhParser.load does not support file-extension .%s; consider using loadx', fext);
        end
        function this = loadx(fn, ext, varargin)
            ip = inputParser;
            addParameter(ip, 'N', mlpet.Resources.instance.defaultN, @islogical);
            parse(ip, varargin{:});
            
            if (~lstrfind(fn, ext))
                if (~strcmp('.', ext(1)))
                    ext = ['.' ext];
                end
                fn = [fn ext];
            end
            assert(lexist(fn, 'file'));
            [pth, fp, fext] = filepartsx(fn, ext); 
            this = mlfourdfp.IfhParser.loadText(fn);
            this.filepath_   = pth;
            this.fileprefix_ = fp;
            this.filesuffix_ = fext;
            this.N_ = ip.Results.N;
        end
        function s = strrep4regexp(s)
            s = strrep(s,'[','\[');
            s = strrep(s,']','\]');
            s = strrep(s,'(','\(');
            s = strrep(s,')','\)');
        end
    end
    
	methods
        function s = asstruct(this, varargin)
            %% ASSTRUCT
            %  @return struct s with fields:
            %  version_of_keys
            %  number_format
            %  conversion_program
            %  name_of_data_file
            %  number_of_bytes_per_pixel
            %  imagedata_byte_order
            %  orientation
            %  number_of_dimensions
            %  matrix_size
            %  scaling_factor
            %  mmppix, if ~this.N
            %  center, if ~this.N
            %  @return s.version_of_keys >= this.MIN_VERSION.
            
            if (~isempty(this.denovo_))
                s = this.denovo_;
                return
            end
            
            assert(this.rightSideNumeric('version of keys') >= this.MIN_VERSION);            
            s = struct( ...
                'version_of_keys', this.rightSideNumeric('version of keys'), ...
                'number_format', this.rightSideChar('number format'), ...
                'conversion_program', this.rightSideChar('conversion program'), ...
                'name_of_data_file', this.nameOfDataFile, ...
                'patient_ID', this.rightSideChar('patient ID'), ...
                'date', this.rightSideChar('date'), ...
                'number_of_bytes_per_pixel', this.rightSideNumeric('number of bytes per pixel'), ...
                'imagedata_byte_order', this.rightSideChar('imagedata byte order'), ...
                'orientation', this.rightSideNumeric('orientation'), ...
                'number_of_dimensions', this.rightSideNumeric('number of dimensions'), ...
                'matrix_size', this.matrixSize, ...
                'global_minimum', this.rightSideNumeric('global minimum'), ...
                'global_maximum', this.rightSideNumeric('global maximum'), ...
                'scaling_factor', this.scalingFactor, ...
                'slice_thickness', this.sliceThickness);
            if (~this.N_), ...
                s.mmppix = this.mmppix;
                s.center = this.center;     
            end
        end
        function n = matrixSize(this)
            idx = 1;
            n = [];
            while (idx < length(this.cellContents))
                [rsn,idx] = this.rightSideNumeric('matrix size [', idx);
                if (isempty(rsn))
                    return
                end
                idx = idx + 1;
                n = [n rsn]; %#ok<AGROW>
            end
        end
        function s = nameOfDataFile(this)
            [~,s] = myfileparts(this.rightSideChar('name of data file'));
        end
        function     save(this, varargin)
            
            ip = inputParser;
            addOptional(ip, 'client', this);
            parse(ip, varargin{:});
            
            fid = fopen(this.fqfilename, 'w'); % overwrite
            fprintf(fid, 'INTERFILE\t:=\n');
            keys = this.SUPPORTED_KEYS;
            str = this.asstruct;
            for ik = 1:length(keys)
                
                if (strcmp(keys{ik}, 'conversion_program'))
                    fprintf(fid, 'conversion program\t:= %s\n', class(ip.Results.client));
                    continue
                end
                if (strcmp(keys{ik}, 'name_of_data_file'))
                    fprintf(fid, 'name of data file\t:= %s\n', this.fileprefix);
                    continue
                end
                if (strcmp(keys{ik}, 'mmppix'))
                    if (~isempty(this.mmppix) && ~this.N_)
                        m = str.mmppix;
                        fprintf(fid, 'mmppix\t:=\t%9.6f %9.6f %9.6f\n', m(1), m(2), m(3));
                    end
                    continue
                end  
                if (strcmp(keys{ik}, 'center')) 
                    if(~isempty(this.center) && ~this.N_)
                        c = str.center;
                        fprintf(fid, 'center\t:=\t%9.4f %9.4f %9.4f\n', c(1), c(2), c(3));
                    end
                    continue
                end   
                
                if (ischar(str.(keys{ik})))
                    fprintf(fid, '%s\t:= %s\n', strrep(keys{ik}, '_',' '), str.(keys{ik}));
                    continue
                end
                if (isnumeric(str.(keys{ik})) && length(str.(keys{ik})) == 1)
                    fprintf(fid, '%s\t:= %g\n', strrep(keys{ik}, '_',' '), str.(keys{ik}));
                    continue
                end
                if (strcmp(keys{ik}, 'scaling_factor'))
                    this.fprintfMulti_f(fid, [strrep(keys{ik}, '_',' ') ' (mm/pixel)'], str.(keys{ik}));
                    continue
                end   
                if (isnumeric(str.(keys{ik})) && length(str.(keys{ik})) > 1)
                    this.fprintfMulti_g(fid, strrep(keys{ik}, '_',' '), str.(keys{ik}));
                    continue
                end
                error('mlfourdfp:guardingIfsFailed', 'in IfhParser.save');
            end
            fprintf('\n');
            fclose(fid);
        end
        function n = scalingFactor(this)
            idx = 1;
            n = [];
            while (idx <= length(this.cellContents))
                [rsn,idx] = this.rightSideNumeric('scaling factor (mm/pixel) [', idx);
                if (isempty(rsn))
                    return
                end
                idx = idx + 1;
                n = [n rsn]; %#ok<AGROW>
            end
        end
        function n = sliceThickness(this)
            n = this.rightSideNumeric('slice thickness (mm/pixel)');
        end
        function n = mmppix(this)
            n = this.rightSideNumerics('mmppix');
        end
        function n = center(this)
            n = this.rightSideNumerics('center');
        end
        function [contnt,idx] = findNextCell(this, fieldName, idx0)
            %  @return idx | cellContents{idx} contains the fieldName.
            
            assert(ischar(fieldName));
            assert(isnumeric(idx0));
            contnt = [];
            idx = idx0;
            for c = idx0:length(this.cellContents_) 
                if (lstrfind(this.cellContents_{c}, fieldName))
                    contnt = this.cellContents_{c};
                    idx = c;
                    break
                end
            end
            %if (isempty(contnt))
            %    error('mlio:endOfFile', 'IfhParser.findNextCell found nothing more'); end
        end	
        function [contntCells,idx] = findNextNCells(this, fieldName, idx0, N)
            assert(ischar(fieldName));
            assert(isnumeric(idx0));
            assert(isnumeric(N));
            contntCells = cell(1,N);
            idx = idx0;
            for c = idx0:length(this.cellContents_) 
                if (lstrfind(this.cellContents_{c}, fieldName))
                    for n = 1:N
                        contntCells{n} = this.cellContents_{c+n-1};
                    end
                    idx = c+N-1;
                    break
                end
            end
            if (isempty(contntCells))
                error('mlio:endOfFile', 'IfhParser.findNextNCells found nothing more'); end
        end	
        function [ch,idx1] = rightSideChar(this, fieldName, varargin)
            p = inputParser;
            addRequired(p, 'fieldName', @ischar);
            addOptional(p, 'idx0', 1,   @isnumeric);
            parse(p, fieldName, varargin{:});
            
            [line,idx1] = this.findNextCell(fieldName, p.Results.idx0);
            if (isempty(line))
                ch = ''; 
                return
            end
            names = regexp(line, sprintf('%s\\s+:=\\s*(?<value1>\\S+)', fieldName), 'names');
            ch    = strtrim(names.value1);
        end
        function [nv,idx1] = rightSideNumeric(this, fieldName, varargin)
            p = inputParser;
            addRequired(p, 'fieldName', @ischar);
            addOptional(p, 'idx0', 1,   @isnumeric);
            parse(p, fieldName, varargin{:});
            
            [line,idx1] = this.findNextCell(fieldName, p.Results.idx0);
            if (isempty(line))
                nv = []; 
                return
            end
            names = regexp(line, sprintf('%s\\S*\\s+:=\\s*(?<value1>%s)', ...
                    this.strrep4regexp(fieldName), this.ENG_PATT_LOW), 'names');
            nv    = str2double(strtrim(names.value1)); 
        end
        function [nv,idx1] = rightSideNumerics(this, fieldName, varargin)
            p = inputParser;
            addRequired(p, 'fieldName', @ischar);
            addOptional(p, 'idx0', 1,   @isnumeric);
            parse(p, fieldName, varargin{:});
            
            [line,idx1] = this.findNextCell(fieldName, p.Results.idx0);
            if (isempty(line))
                nv = []; 
                return
            end
            names = regexp(line, sprintf('%s\\S*\\s+:=\\s*(?<value1>\\S+\\s+\\S+\\s+\\S+)', ...
                    this.strrep4regexp(fieldName)), 'names');
            nv    = str2num(strtrim(names.value1)); %#ok<ST2NM>
        end
        
        function this = IfhParser(varargin)            
            ip = inputParser;
            addParameter(ip, 'N', mlpet.Resources.instance.defaultN, @islogical);
            parse(ip, varargin{:});
            
            this.filesuffix = this.IFH_EXT;
            this.N_ = ip.Results.N;
        end
    end 
    
    %% PROTECTED
    
    properties (Access = 'protected')
        denovo_
        N_ 
    end
    
    methods (Static, Access = 'protected')
        function        fprintfMulti_f(fid, key, val)
            assert(isnumeric(val));
            for iv = 1:length(val)
                fprintf(fid, '%s [%i]\t:= %8.6f\n', key, iv, val(iv));
            end
        end
        function        fprintfMulti_g(fid, key, val)
            assert(isnumeric(val));
            for iv = 1:length(val)
                fprintf(fid, '%s [%i]\t:= %g\n', key, iv, val(iv));
            end
        end
        function this = loadText(fn)
            import mlfourdfp.*;
            this = IfhParser;
            this.cellContents_ = IfhParser.textfileToCell(fn);
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

