classdef ImgRecParser < handle & mlio.AbstractParser
	%% IMGRECPARSER parses numerical values to the right or left of a text field-name.
    %  For more fine-grained parsing features, see TextParser.

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.5.0.197613 (R2015a) 
 	%  $Id$   
    
    properties (Constant)
        IMGREC_EXT = {'.4dfp.img.rec' };
    end
    
	methods (Static)
        function this = load(fn)
            assert(lexist(fn, 'file'));
            [pth, fp, fext] = myfileparts(fn); 
            if (lstrfind(fext, mlfourdfp.ImgRecParser.IMGREC_EXT) || ...
                isempty(fext))
                this = mlfourdfp.ImgRecParser.loadText(fn); 
                this.filepath_   = pth;
                this.fileprefix_ = fp;
                this.filesuffix_ = fext;
                return 
            end
            error('mlio:unsupportedParam', ...
                'ImgRecParser.load does not support file-extension .%s; consider using loadx', fext);
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
            this = mlfourdfp.ImgRecParser.loadText(fn);
            this.filepath_   = pth;
            this.fileprefix_ = fp;
            this.filesuffix_ = fext;
        end
    end
    
	methods
        function srow = commonSform(this)
            [~,line] = this.findFirstCell('common sform[0]:');
            if (isempty(line))
                srow = [];
                return
            end
            posColon = strfind(this.cellContents_{line}, ':');
            content  = this.cellContents_(line:line+2);
            content  = cellfun(@(x) x(posColon+1:end), content, 'UniformOutput', false);
            srow     = [str2num(content{1}); str2num(content{2}); str2num(content{3})]; %#ok<ST2NM>
            srow(1,:) = -srow(1,:);
        end
        function cntnt = extractLinesByRegexp(this, re)
            idxs = regexp(this.cellContents_, re);
            cntnt = {};
            for c = 1:length(this.cellContents_)
                if (~isempty(idxs{c}))
                    cntnt = [cntnt this.cellContents_{c}]; %#ok<AGROW>
                end
            end
        end
        function [contnt,idx] = findNextCell(this, fieldName, idx0)
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
            if (isempty(contnt))
                error('mlio:endOfFile', 'ImgRecParser.findNextCell found nothing more'); end
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
                error('mlio:endOfFile', 'ImgRecParser.findNextNCells found nothing more'); end
        end	
    end 
    
    %% PROTECTED
    
    methods (Static, Access = 'protected')
        function this = loadText(fn)
            import mlfourdfp.*;
            this = ImgRecParser;
            this.cellContents_ = ImgRecParser.textfileToCell(fn);
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

