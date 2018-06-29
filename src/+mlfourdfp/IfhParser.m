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
        IFH_EXT = {'.ifh' };
        MIN_VERSION = 3.3
    end
    
	methods (Static)
        function this = load(fn)
            assert(lexist(fn, 'file'));
            [pth, fp, fext] = myfileparts(fn); 
            if (lstrfind(fext, mlfourdfp.IfhParser.IFH_EXT) || ...
                isempty(fext))
                this = mlfourdfp.IfhParser.loadText(fn); 
                this.filepath_   = pth;
                this.fileprefix_ = fp;
                this.filesuffix_ = fext;
                return 
            end
            error('mlio:unsupportedParam', ...
                'IfhParser.load does not support file-extension .%s; consider using loadx', fext);
        end
        function this = loadx(fn, ext)
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
        end
        function s = strrep4regexp(s)
            s = strrep(s,'[','\[');
            s = strrep(s,']','\]');
            s = strrep(s,'(','\(');
            s = strrep(s,')','\)');
        end
    end
    
	methods
        function s = asstruct(this)
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
            %  mmppix
            %  center
            %  @return s.version_of_keys >= this.MIN_VERSION.
            
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
                'slice_thickness', this.sliceThickness, ...
                'mmppix', this.mmppix, ...
                'center', this.center);     
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
    end 
    
    %% PROTECTED
    
    methods (Static, Access = 'protected')
        function this = loadText(fn)
            import mlfourdfp.*;
            this = IfhParser;
            this.cellContents_ = IfhParser.textfileToCell(fn);
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

