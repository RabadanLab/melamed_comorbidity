function numMat = cell2num(cellArray)
  numMat = zeros(size(cellArray));
  for i = 1:size(cellArray,2)
    numMat(:,i) = str2double(cellArray(:,i));
  end