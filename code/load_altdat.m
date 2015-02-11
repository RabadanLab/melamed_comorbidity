function [altdat, sel, rr] = load_altdat(dirname) 
%alts = {'amp_peak','del_peak','mutation','peak_mut'};
alts = {'peak_mut'};
altdat = struct();

res = parseText([dirname '/' alts{1} '.xls'],'colname',true,'nlines',3);
CUT = find(ismember(res.colname,'relative_risk')) - 1;

for alt = 1:length(alts)
  res = parseText([dirname '/' alts{alt} '.xls'],'colname',true);

  %res.footer = cell2num(res.text(end,4:end));
  thedat = cell2num(res.text(:,(CUT + 1):end));
  thetext = res.text(:,1:CUT);
  res.gcolname = res.colname(1:CUT);
  for j = 1:length(res.gcolname)
    res.(res.gcolname{j}) = thetext(:,j);
  end

  res.colname = res.colname((CUT+1):end);
  for j = 1:length(res.colname)
    res.(res.colname{j}) = thedat(:,j);
  end
  res.text = 1;
  altdat.(alts{alt}) = res;

end

sel = altdat.(alts{1}).relative_risk > -1 & altdat.(alts{1}).gene_enrichment > -1;
rr = altdat.(alts{alt}).relative_risk(sel) > 1;
