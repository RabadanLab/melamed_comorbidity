function pdf(ptitle)
  set(gcf,'PaperPositionMode','auto')
  print('-dpdf','-r100',[ptitle '.pdf'])