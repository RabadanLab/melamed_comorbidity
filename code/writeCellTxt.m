function writeCellTxt(strings, fname, delim, colNames, mode_in)
del = '\t'; if nargin>2 del = delim; end
mode = 'wt';
if nargin > 4 mode = mode_in; end
dfo = fopen(fname,mode);
if ~iscellstr(strings(1,:))
  for i = 1:size(strings,2)
    if ~iscellstr(strings(:,i)) % assume everything else is numeric
      strings(:,i) = strtrim(cellstr(num2str([strings{:,i}]')));
    end
  end
end
if nargin > 3 
  % this might take too much memory: strings = [colNames; strings]; end
  fprintf(dfo,['%s' del], colNames{1:(end-1)});
  fprintf(dfo,'%s\n',colNames{end})
end
for i = 1:size(strings,1)
  fprintf(dfo,['%s' del], strings{i,1:(end-1)});
  fprintf(dfo,'%s\n', strings{i,end});
end
fclose(dfo);
