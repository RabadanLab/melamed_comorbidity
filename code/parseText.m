function res = parseText(fn, varargin)
    %fn: file name
    %varargin:
    %   'skip': number of lines to skip at the beginning of the file
    %       (comments etc)
    %   'colname': {true, false}; def=false;
    %   'rowname': {true, false}; def=false; if true, rownames are extracted
    %       from first column
    %   'ncol': number of columns if known
    %   'delimiter': default=tab
    %   'numeric': if the text body is numeric (except for the rowname, if
    %       exists); default=false;
    %   'nlines':  how many lines to read
    %
    %Note: to make sure it reads file correctly, put NaN for ''
    %
    %return:
    %   res.rowname, if rowname=true
    %   res.colname, if colname=true
    %   res.text
    %
    
    %default
    delimiter = 9; %tab
    colname = false;
    rowname = false;
    skip = 0;
    ncol = 0;
    numeric = false;
    nlines = inf;
    
    i = 1;
    while i < length(varargin)
      switch lower(varargin{i})
        case {'colname','rowname','numeric'}
          eval(sprintf('%s = logical(%d);',varargin{i},varargin{i+1}));
        case {'skip','ncol'}
          eval(sprintf('%s = %d;',varargin{i},varargin{i+1}));
        case {'delimiter'}
          if strcmp(varargin{i+1},'\t')
            delimiter = 9;
          else
            delimiter = varargin{i+1};
          end
        case {'nlines'}
          nlines = varargin{i+1};;
        otherwise
          error('Unknown option %s.\n',varargin{i});
      end
      i = i + 2;
    end
    
    fid = fopen(fn,'rt');    
    %skip lines if any
    i = 0;
    while i < skip 
        line = fgetl(fid);
        i = i + 1;
    end
    if colname
        hline = fgetl(fid);
        res.colname = {};
        [res.colname{end+1} t] = strtok(hline,delimiter);
        while length(t) > 0
            [res.colname{end+1} t] = strtok(t, delimiter);
        end
        res.colname = res.colname';
        ncol = size(res.colname,1);
        if numeric ncol = ncol + 1; end
    end
    
    if numeric
        formattail = ' %f';     res.text = [];
    else
        formattail = ' %s';     res.text = {};
    end
    if rowname || ~numeric %at least the file has one column
        format = '%s';
    else
        format = '%f';
    end
    if ncol == 0
        fpos = ftell(fid);
        line = fgetl(fid);
        offset = fpos - ftell(fid);        
        fseek(fid, offset, 'cof');
        %decide the number of columns        
        [a t] = strtok(line, delimiter);
        if length(a) > 0            
            ncol = ncol + 1;
        end        
        while length(t) > 0
            [a t] = strtok(t, delimiter);
            format = [format formattail];
            ncol = ncol + 1;
        end
    else        
        for i = 2:ncol
            format = [format formattail];
        end
    end
    if isscalar(delimiter), delimiter = '\t'; end
    text = [];
    if nlines < Inf
      text = textscan(fid, format, nlines, 'delimiter',delimiter);    
    else
      text = textscan(fid, format, 'delimiter',delimiter,'bufsize',49140);    
    end
    fclose(fid);
    
    if rowname
        res.rowname = text{1};
        text = text(2:end);
    end
    
    for i = 1:length(text)
        res.text(:,i) = text{i};
    end

    
    
    