def draw_rect(dic, name, ax, xRange, yRange, colormap,
              normalize=None, logrange=1e10,
              edge_color=None, edge_width=None, draw_box_only=False,
              rasterized=False):
  from  matplotlib.collections import PolyCollection
  import numpy as np
  poly_collec = []
  facecolors = []
  edgecolors = []
  linewidths = []
  nlen = dic['rmin'].size

  maxval = np.nanmax(dic[name])
  minval = np.nanmin(dic[name])

  if normalize == 'linear':
    def norm_f(v):
      return (v - minval) / (maxval - minval)
  elif normalize == 'log':
    logmaxval = np.log(maxval)
    if minval * logrange >= maxval:
      logminval = np.log(minval)
    else:
      logminval = np.log(maxval / logrange)
    def norm_f(v):
      if v <= 0:
        return 0.0
      return (np.log(v) - logminval) / (logmaxval - logminval)
  else:
    norm_f = normalize

  for i in range(nlen):
      # Make rectangle
      x1 = dic['rmin'][i]
      y1 = dic['zmin'][i]
      x2 = dic['rmax'][i]
      y2 = dic['zmax'][i]
      if x1 > xRange[1] or y1 > yRange[1]:
          continue
      if x2 < xRange[0] or y2 < yRange[0]:
          continue
      pxy = np.zeros((5,2))
      pxy[0, :] = [x1, y1]
      pxy[1, :] = [x2, y1]
      pxy[2, :] = [x2, y2]
      pxy[3, :] = [x1, y2]
      pxy[4, :] = [x1, y1]
      #
      # Calculate the color
      thiscolor = colormap(norm_f(dic[name][i]))
      #
      # Draw the rectangle filled with color
      if draw_box_only:
        facecolor = 'none'
        edgecolor = 'black'
        linewidth = 0.1
      else:
        facecolor = thiscolor
        if edge_color != None or edge_width != None:
            edgecolor = edge_color if edge_color != None else 'black'
            edgewidth = edge_width if edge_width != None else 0.1
        else:
            edgecolor = thiscolor
            edgewidth = 0.00001
      poly_collec.append(pxy)
      facecolors.append(facecolor)
      edgecolors.append(edgecolor)
      linewidths.append(edgewidth)

  ax.add_collection(PolyCollection(poly_collec, edgecolors=edgecolors,
      facecolors=facecolors, linewidths=linewidths, rasterized=rasterized))
  
