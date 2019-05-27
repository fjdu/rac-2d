def plot_a_pannel(isub, dic, idx, name, x1, x2, xr, yr, minval=1e-50,
                  xscale='linear', yscale='log', linewidth=2):
    import matplotlib.pyplot as plt
    plt.subplot(*isub); isub[2] += 1
    y = dic[name][idx]
    stepplot(x1, x2, y, minval=minval, linewidth=linewidth)
    plt.gca().set_xlim(xr)
    plt.gca().set_ylim(yr)
    plt.gca().set_xscale(xscale)
    plt.gca().set_yscale(yscale)
    set_axis_format(plt.gca(), xscale=xscale, yscale=yscale)
    plt.text(0.5, 0.97, name, transform = plt.gca().transAxes, fontsize=15,
         horizontalalignment='center', verticalalignment='top')


class query_state:
    def __init__(self):
        self.saved = {}

    def query_item(self, d, rho, z):
        found = False
        nr = len(d['rmin'])
        for i in range(nr):
            if (d['rmin'][i] <= rho <= d['rmax'][i]) and (d['zmin'][i] <= z <= d['zmax'][i]):
                found = True
                break
        if found:
            return i
        else:
            return None

    def query(self, d, rho, z, key):
        if (rho, z) not in self.saved:
            idx = self.query_item(d, rho, z)
            self.saved[(rho, z)] = idx
        else:
            idx = self.saved[(rho, z)]
        if idx != None:
            return d[key][idx]
        else:
            return 0.0


def write_spherical_coordinates(filename, r_grid, theta_grid, phi_grid):
    with open(filename, 'w') as f:
        f.write()
        f.write()
    


def set_axis_format(ax, xscale=None, yscale=None, labelfontsize=25, tickfontsize=20, majorgridon=True, minorgridon=True,
                    onlyfont=False, graygrid=False, xmingrid=False, xmajgrid=False, ymingrid=False, ymajgrid=False,
                    color_major=(0.7,0.8,0.7), color_minor=(0.5,0.9,0.7)):
    import numpy as np
    from matplotlib.ticker import AutoMinorLocator, LogLocator
    if xscale == None:
      xscale = ax.get_xscale()
    if yscale == None:
      yscale = ax.get_yscale()
    if not onlyfont:
      if xscale == 'linear':
          ax.xaxis.set_minor_locator(AutoMinorLocator(10))
      if xscale == 'log':
          xt = ax.get_xticks()
          ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.linspace(1.0, 10.0**np.log10(xt[1]/xt[0]), 10)))
      if yscale == 'linear':
          ax.yaxis.set_minor_locator(AutoMinorLocator(10))
      if yscale == 'log':
          yt = ax.get_yticks()
          ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.linspace(1.0, 10.0**np.log10(yt[1]/yt[0]), 10)))
      if graygrid:
          color_major=(0.5,0.5,0.5)
          color_minor=(0.8,0.8,0.8)
      if majorgridon:
          xmajgrid = True
          ymajgrid = True
      if minorgridon:
          xmingrid = True
          ymingrid = True
      if xmajgrid:
          ax.xaxis.grid(which='major', color=color_major, linewidth=1,   linestyle='-')
      if xmingrid:
          ax.xaxis.grid(which='minor', color=color_minor, linewidth=0.5,   linestyle='-')
      if ymajgrid:
          ax.yaxis.grid(which='major', color=color_major, linewidth=1, linestyle='-')
      if ymingrid:
          ax.yaxis.grid(which='minor', color=color_minor, linewidth=0.5, linestyle='-')
      #
      ax.set_axisbelow(True)
    #
    ax.xaxis.label.set_fontsize(labelfontsize)
    ax.yaxis.label.set_fontsize(labelfontsize)
    for label in ax.get_xticklabels():
        label.set_fontsize(tickfontsize)
    for label in ax.get_yticklabels():
        label.set_fontsize(tickfontsize)



def plot_model_results_map(dic, lis, panel_pos=(0.10,0.15,0.75,0.8), bgcolor='none',
                           hcmode=False, hc_list=None, hc_colors=None, useSecond=False,
                           lgd_alpha=0.5, lgd_bbox=None, showlgd=True,
                           draw_rectangles = True, fig_prev=None, ax_prev=None,
                           ax_cbar=None,cbar_orien=u'vertical', cbar_tick_loc='auto',
                           xRange = None, yRange = None,
                           xscale='linear', yscale='linear',
                           xtitle='r (AU)', ytitle='z (AU)',
                           showColorbar=True, rasterized=True,
                           forcegte=False, forcelte=False,
                           figsize = None, display_dic=None, show_name=True, name_position=(0.5,0.98),
                           labelfontsize=20, tickfontsize=16, textfontsize=20,
                           majorgridon=True, minorgridon=True,
                           xmingrid=False, xmajgrid=False, ymingrid=False, ymajgrid=False,
                           cbarlabel='', useoldstylecbar=False,
                           linthreshy=1e-3, graygrid=False,
                           all_in_one_fig=False, pansepxfrac=0.2, pansepyfrac=0.09,
                           edge_color=None,
                           edge_width=None,
                           draw_ray_to=None,
                           draw_ray_slope_range=None,
                           draw_ray_r_range=None,
                           draw_box_only=False,
                           nx = 100, ny = 100, nlev = 100):
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import AutoLocator, LogLocator
    import os
    #import my_script as my
    from draw import my_script as my
    from scipy.interpolate import griddata
    from  matplotlib.collections import PolyCollection
    #
    polygonzorder = 0
    #
    minr = np.min(dic['rmin'])
    maxr = np.max(dic['rmax'])
    minz = np.min(dic['zmin'])
    maxz = np.max(dic['zmax'])
    if xRange == None:
        xRange = (minr, maxr)
    if yRange == None:
        if yscale == 'log':
          minz = np.min(dic['zmax'])*0.8
          yRange = (minz, maxz)
        else:
          yRange = (minz, maxz)
    #
    if figsize == None:
        figsize = (10, 8)
    #
    if hcmode:
      filelen = dic['rmin'].shape[0]
      hc_prep_hc(dic, hc_list)
    else:
      filelen = dic[lis[0][0]['name']].shape[0]
    #
    matplotlib.rcParams['axes.linewidth'] = 1
    rcpara_tex = matplotlib.rcParams['text.usetex']
    matplotlib.rcParams['text.usetex'] = True
    #
    if hcmode:
      hc_used = [False for j in range(len(hc_list))]
    #
    if all_in_one_fig:
      fig = plt.figure(figsize=figsize)
      npanx = np.ceil(np.sqrt(len(lis)))
      npany = np.ceil(float(len(lis))/npanx)
      xmarginleft = 0.1
      ymarginlower = 0.1
      panwx = (0.99 - xmarginleft)  / npanx
      panwy = (0.99 - ymarginlower) / npany
      pansepx = panwx * pansepxfrac
      pansepy = panwy * pansepyfrac
      panwx -= pansepx
      panwy -= pansepy
    #
    ax_s = {}
    #
    igroup = 0
    for item_group in lis:
        #
        igroup += 1
        #
        # All the items in a group in the list are plotted in the same pannel
        if not all_in_one_fig:
          if fig_prev == None:
            fig = plt.figure(figsize=figsize)
          else:
            fig = fig_prev
        #
        len_group = len(item_group)
        iitem = 0
        for item in item_group:
            name = item['name']
            if 'linearscale' in item:
                use_linear_scale = item['linearscale']
            else:
                use_linear_scale = False
            #
            if 'cmap' in item:
                colormap = item['cmap']
            else:
                colormap = cm.rainbow
            #
            if not hcmode and 'cdict' not in item:
                maxval = np.nanmax(dic[name])
                minval = np.nanmin(dic[name])
                idx_nonzero = np.where(dic[name] > 0.0)
                if len(idx_nonzero[0]) > 0:
                  minval_nonzero = np.nanmin(dic[name][idx_nonzero])
                else:
                  minval_nonzero = None
                print('{:12s}{:s}{:12.4e}{:12.4e}'.format(name, ' Max, min ', maxval, minval))
            min_omit = False
            max_omit = False
            if 'vr' in item:
                if len(item['vr']) >= 2:
                    if minval < item['vr'][0]:
                      min_omit = True
                    if maxval > item['vr'][1]:
                      max_omit = True
                    minval = item['vr'][0]
                    minval_nonzero = item['vr'][0]
                    maxval = item['vr'][1]
            if not hcmode and 'cdict' not in item:
              if not use_linear_scale:
                log_max = np.log10(maxval)
                log_min = np.log10(max(minval, minval_nonzero, maxval/1e20))
                minval = 10**(log_min)
            #
            iitem += 1
            #
            if all_in_one_fig:
                kx = np.mod(igroup-1, npanx)
                ky = npany - np.ceil(igroup/npanx)
                xleft  = xmarginleft  + panwx * kx + pansepx * kx
                ylower = ymarginlower + panwy * ky + pansepy * ky
                pos = (xleft, ylower, panwx, panwy)
            else:
                pos = panel_pos
            if (ax_prev == None) and (iitem == 1):
                ax = fig.add_axes(pos,
                              xlabel=xtitle if iitem == len_group else '',
                              ylabel=ytitle, facecolor=bgcolor,
                              #rasterized=rasterized,
                              autoscalex_on=False, autoscaley_on=False,
                              xscale=xscale, yscale=yscale,
                              xlim=xRange, ylim=yRange)
                ax_s.update({name: ax})
                if yscale == 'symlog':
                    ax.set_yscale('symlog', linthreshy=linthreshy)
                if all_in_one_fig:
                    if kx > 0:
                        ax.set_ylabel('')
                        ax.set_yticklabels([])
                    if ky > 0:
                        ax.set_xlabel('')
                        ax.set_xticklabels([])
            else:
                ax = ax_prev
            if not ax:
              ax = plt.gca()
            #
            if show_name:
              if 'name_disp' in item:
                  name_disp = item['name_disp']
              else:
                  name_disp = name
                  name_disp = name_disp.replace('_', '\\_')
                  if display_dic != None:
                    if name in display_dic:
                      name_disp = display_dic[name]
            else:
              name_disp = ''
            plt.text(name_position[0], name_position[1], name_disp, transform = plt.gca().transAxes,
                 fontsize=textfontsize, horizontalalignment='center', verticalalignment='top')
            set_axis_format(ax, xscale=xscale, yscale=yscale, graygrid=graygrid,
                labelfontsize=labelfontsize, tickfontsize=tickfontsize,
                majorgridon=majorgridon, minorgridon=minorgridon,
                xmingrid=xmingrid, xmajgrid=xmajgrid, ymingrid=ymingrid, ymajgrid=ymajgrid)
            ax.set_axisbelow(False)
            #
            if draw_rectangles or hcmode:
                poly_collec = []
                facecolors = []
                edgecolors = []
                linewidths = []
                for i in range(filelen):
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
                    if hcmode:
                        hc_vals = np.array([dic[list(hc_list[j].keys())[0]][i] for j in range(len(hc_list))])
                        if useSecond:
                          imax = hc_vals.argsort()[-2]
                        else:
                          imax = hc_vals.argmax()
                        thiscolor = hc_colors[imax]
                        hc_used[imax] = True
                    elif 'cdict' in item:
                        thiscolor = item['cdict'][dic[name][i]]
                    else:
                        val = dic[name][i]
                        if use_linear_scale:
                            sca_col = (val - minval) / (maxval - minval)
                        else:
                            if val <= 0.0:
                                sca_col = np.nan
                            else:
                                sca_col = (np.log10(val) - log_min) / (log_max-log_min)
                        thiscolor = colormap(sca_col)
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
                    #
                    if draw_ray_to != None:
                        draw_ray = True
                        if draw_ray_slope_range != None:
                          slope = (y1+y2) / (x1+x2)
                          if slope < draw_ray_slope_range[0] or slope > draw_ray_slope_range[1]:
                            draw_ray = False
                        if draw_ray_r_range != None and draw_ray:
                          r = sqrt((y1+y2)**2 + (x1+x2)**2)*0.5
                          if r < draw_ray_r_range[0] or r > draw_ray_r_range[1]:
                            draw_ray = False
                        if draw_ray:
                          ax.plot([0.5*(x1+x2), draw_ray_to[0]],
                                  [0.5*(y1+y2), draw_ray_to[1]], linewidth=0.2, color=(1.0,0.7,0.7))
                          ax.plot([0.5*(x1+x2)],
                                  [0.5*(y1+y2)], marker='.', markersize=1, color=(1.0,1.0,1.0))
                ##########
                ax.add_collection(PolyCollection(poly_collec, edgecolors=edgecolors,
                    facecolors=facecolors, linewidths=linewidths, rasterized=rasterized))
            else:
                x = dic['rmax']
                y = dic['zmin'] #0.5 * (dic['zmin'] + dic['zmax'])
                xi = np.linspace(xRange[0], xRange[1], nx)
                yi = np.linspace(yRange[0], yRange[1], ny)
                xi, yi = np.meshgrid(xi, yi)
                z = np.zeros_like(dic[name]) - 100.0
                if len(idx_nonzero) > 0:
                  z[idx_nonzero] = np.log10(dic[name][idx_nonzero])
                else:
                  return
                zi = griddata((x, y), z, (xi, yi), method='cubic')
                zi[zi < log_min] = log_min
                zi[zi > log_max] = log_max
                ax.contour(xi, yi, zi, nlev,
                           cmap=colormap, linewidth=0.001,
                           vmin=log_min, vmax=log_max)
                ax.contourf(xi, yi, zi, nlev, cmap=colormap, vmin=log_min, vmax=log_max)
                if 'levels' in item:
                    CS = ax.contour(xi, yi, zi, levels=np.log10(item['levels']),
                                    linestyles='solid', linewidth=1,
                                    colors='white', vmin=log_min, vmax=log_max)
                    plt.clabel(CS, inline=1, fontsize=10)                    
            #
            if rasterized == True:
                ax.set_rasterization_zorder(polygonzorder+1)
            #
            if 'legendFunc' in item:
                item['legendFunc'](ax)
            #
            if showColorbar and (not hcmode) and (not draw_box_only) and \
               (not ('showColorbar' in item and item['showColorbar'] == False)):
                ntick_cbar = 5
                #ax1 = plt.subplot2grid((len_group, 15), (iitem-1, 13), colspan=1)
                pwidth = pos[2]
                pright = pos[0] + pwidth
                cbarwidth = pwidth * 0.03
                if ax_cbar == None:
                  pos = (pright+cbarwidth*0.3, pos[1], cbarwidth, pos[3])
                  ax1 = fig.add_axes(pos)
                else:
                  ax1 = ax_cbar
                if use_linear_scale:
                    if useoldstylecbar:
                        tickvals = np.linspace(minval, maxval, num=ntick_cbar)
                    ax1.set_xscale('linear')
                    ax1.set_xlim((minval, maxval))
                else:
                    if useoldstylecbar:
                        tickvals = np.logspace(log_min, log_max, num=ntick_cbar)
                    ax1.set_xscale('log')
                    ax1.set_xlim((10.0**log_min, 10.0**log_max))
                if useoldstylecbar:
                    ticks = np.linspace(0.0, 1.0, ntick_cbar)
                else:
                    tickvals = ax1.xaxis.get_ticklocs()
                    ticks = [ax1.transAxes.inverted().transform(\
                            ax1.transData.transform(np.array([(tic,0.1)])))[0][0] for tic in tickvals]
                ticklabels = [my.num2str4fig(tickval) for tickval in tickvals]
                ticks, ticklabels = zip(*filter(lambda x: 0<=x[0]<=1, zip(ticks, ticklabels)))
                if len(ticks) < 5:
                    if use_linear_scale:
                      L = AutoLocator()
                    else:
                      L = LogLocator(base=10, subs=np.linspace(1.0,10.0,10))
                    tkvals = L.tick_values(minval, maxval)
                    minticks = [ax1.transAxes.inverted().transform(\
                            ax1.transData.transform(np.array([(tic,0.1)])))[0][0] for tic in tkvals]
                    minticks = [_ for _ in minticks if 0 <= _ <= 1]
                ax1.set_xscale('linear')
                #if min_omit or forcelte:
                #  ticklabels[0] = '$\leq$' + ticklabels[0]
                #if max_omit or forcegte:
                #  ticklabels[-1] = '$\geq$' + ticklabels[-1]
                cbar = matplotlib.colorbar.ColorbarBase(ax1, ticks=ticks, ticklocation=cbar_tick_loc,
                                                        cmap=colormap,
                                                        #values=None,
                                                        norm=matplotlib.colors.Normalize(vmin=0, vmax=1, clip=False),
                                                        orientation=cbar_orien)
                cbar.solids.set_rasterized(True)
                cbar.set_ticklabels(ticklabels)
                if cbarlabel == None:
                    cbarlabel = name_disp
                if len(ticks) < 5:
                    cbar.ax.yaxis.set_ticks(minticks, minor=True)
                    cbar.ax.yaxis.set_tick_params(which='minor', length=5)
                    cbar.ax.yaxis.set_tick_params(which='major', length=7)
                cbar.ax.tick_params(labelsize=tickfontsize) 
                cbar.set_label(cbarlabel, size=textfontsize)
            if hcmode and showlgd:
                p = []
                L = []
                for j in range(len(hc_list)):
                  if hc_used[j]:
                    p.append(Rectangle((0,0), 1, 1, fc=hc_colors[j], ec='none'))
                    L.append(hc_list[j][list(hc_list[j].keys())[0]])
                plt.legend(p, L, framealpha=lgd_alpha, bbox_to_anchor=lgd_bbox)
    matplotlib.rcParams['text.usetex'] = rcpara_tex
    if len(ax_s) == 0:
        return (fig, ax)
    else:
        return (fig, ax_s)




def hc_prep_hc(c, hclist):
    for hc in hclist:
      hck = list(hc.keys())[0]
      if not hck in c:
        if hck[0] == 'h':
          k = 'c' + hck[1:]
          try:
            c.update({hck: -c[k]})
          except KeyError:
            print('Key error: ', k, ' not in keys.')
            raise
        elif hck[0] == 'c':
          k = 'h' + hck[1:]
          try:
            c.update({hck: -c[k]})
          except KeyError:
            print('Key error: ', k, ' not in keys.')
            raise

def overlap_contour(fig, ax, dic, name, levels, nxy=100, linewidth=1, color='black', fontsize=12,
                    method='nearest', fmt='%1.3f', manual=False, nwidth=1, whichx='rmin', whichy='zmin'):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata
    import scipy.signal
    x = dic[whichx]
    y = dic[whichy]
    xRange = ax.get_xlim()
    yRange = ax.get_ylim()
    nx = nxy
    ny = nxy
    dx = xRange[1] - xRange[0]
    dy = yRange[1] - yRange[0]
    xi = np.linspace(xRange[0], xRange[1], nx)
    yi = np.linspace(yRange[0], yRange[1], ny)
    xi, yi = np.meshgrid(xi, yi)
    z = dic[name]
    zi = griddata((x, y), z, (xi, yi), method=method)
    if nwidth > 0:
        nsmooth = 3*nwidth+1
        xt = np.linspace(-nsmooth, nsmooth, 2*nsmooth+1)
        xt, yt = np.meshgrid(xt, xt)
        kernel = np.exp(-(xt**2 + yt**2)/(2.0*nwidth))
        kernel /= kernel.sum()
        zi = scipy.signal.convolve2d(zi, kernel, mode='same', boundary='symm')
    CS = ax.contour(xi, yi, zi, levels, colors=color,
               linewidth=linewidth)
    plt.clabel(CS, inline=1, colors=color, fontsize=fontsize, fmt=fmt, manual=manual)


def plot_model_results_column(dic, lis, rpos, xRange = None, figsize = None,
                              xscale='linear', yscale='log', whichx='zmin', divideby=None,
                              z_label=None, ylabel='', lw=2,
                              labelname=None, drawlegend=True, drawtitle=False,
                              marker='None',
                              ax=None, fig=None,
                              returnax=False,
                              color=None,
                              do_interpol=False):
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import os
    import my_script as my
    from scipy import interpolate
    
    color_list = ['blue', 'red', 'green', 'magenta',
                  (0.7,0.5,0.3), (0.8,0.85,0.3), (0.2,0.2,0.2), (0.5,0.5,0.5)]
    
    idx = np.where(np.logical_and(dic['rmin'] <= rpos, dic['rmax'] >= rpos))[0].tolist()
    
    if whichx == 'zmin':
        if z_label==None:
          z_label = 'z (AU)'
        if divideby == None:
            idx, z_using = zip(*sorted(zip(idx,
                                       dic['zmin'][idx]),
                                   key = lambda t: t[1]))
        else:
            idx, z_using, z_dvd = zip(*sorted(zip(idx,
                                       dic['zmin'][idx],
                                       dic[divideby][idx]),
                                   key = lambda t: t[1]))
            z_using = np.array(z_using) / np.array(z_dvd)
    elif whichx == 'Av':
        z_label = 'log$_{10}$Av'
        idx, z_using = zip(*sorted(zip(idx,
                                       np.log10(dic['Av_ISM'][idx])),
                                   key = lambda t: t[1]))
    else:
        z_label = whichx
        idx, z_using = zip(*sorted(zip(idx,
                                       np.log10(dic[whichx][idx])),
                                   key = lambda t: t[1]))
    
    idx = list(idx)
    minz = np.min(z_using)
    maxz = np.max(z_using)
    
    if xRange == None:
        xRange = (minz, maxz)
        
    if figsize == None:
        figsize = (10, 8)
        
    matplotlib.rcParams['axes.linewidth'] = 1
    #
    ngroup = len(lis)
    if fig == None:
        fig = plt.figure(figsize=figsize)
    isubplt = 0
    for item_group in lis:
        # All the items in a group in the list are plotted in the same pannel
        if 'vr' in item_group[0]:
            yRange = item_group[0]['vr']
        else:
            val_tmp = dic[item_group[0]['name']][idx]
            vmax = np.max(val_tmp)
            vmin = np.min(val_tmp)
            yRange = (max(vmin/1.5, vmax/1e10), vmax * 4)
        if ax == None:
            ax = plt.subplot2grid((ngroup, 1), (isubplt, 0), colspan=1,
                              xlabel=z_label, ylabel=ylabel,
                              autoscalex_on=False, autoscaley_on=False,
                              xscale=xscale, yscale=yscale,
                              xlim=xRange, ylim=yRange)
        if isubplt == 0 and drawtitle:
          plt.text(0.5, 0.98, 'r = {0:7.2f} AU'.format(rpos), transform = ax.transAxes,
             fontsize=25, horizontalalignment='center', verticalalignment='top')
        isubplt += 1
        #
        icolor = 0
        for item in item_group:
            name = item['name']
            #f = interpolate.interp1d(z_using, dic[name][idx])
            #znew = np.linspace(minz, maxz, 50)
            #ax.plot(znew, f(znew), linestyle='-', label=name,
            #        color=color_list[icolor%len(color_list)], linewidth=5)
            if labelname==None:
                labelname=name
                if 'display_dic' in globals():
                    if name in display_dic:
                        labelname = display_dic[name]
            if color == None:
                color = color_list[icolor%len(color_list)]
            if do_interpol:
                ax.plot(*loglin_interpol(z_using, dic[name][idx], xRange=xRange, num=100, method='linear'),
                        label=labelname,
                        linewidth=lw, linestyle='-', marker=marker, markersize=5,
                        color=color)
            else:
                ax.plot(z_using, dic[name][idx], label=labelname,
                        linewidth=lw, linestyle='-', marker=marker, markersize=5,
                        color=color)
            icolor += 1
        if drawlegend:
            lgd = ax.legend(loc='lower left', bbox_to_anchor=(0.8, 0.15), prop={'size':15},
                        fancybox=False, shadow=False, ncol=1)
            lgd.get_frame().set_alpha(0.5)
    #plt.tight_layout()
    if returnax:
        return (fig, ax)




def plot_model_results_radial(dic, lis, zpos, xRange = None, figsize = None,
        xscale='linear', yscale='log', returndata=True, noLegend=False, noTitle=False,
        color=None, label=None, display_dic=None,
        ax=None, returnaxis=False, set_ax_format=False, linestyle='-', lw=3, markersize=5,
        do_interpol=False, interpolNum=200,
        dosmooth=False, winwidth=None,
        whichx='rmin'):
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import os
    import my_script as my
    from scipy import interpolate
    
    color_list = ['blue', 'red', 'green', 'magenta',
                  (0.7,0.5,0.3), (0.8,0.85,0.3), (0.2,0.2,0.2), (0.5,0.5,0.5)]

    idx = np.where(np.logical_and(dic['zmin'] <= zpos, dic['zmax'] >= zpos))[0].tolist()

    if whichx == 'rmin':
      x_label = 'r (AU)'
      idx, x_using = zip(*sorted(zip(idx,
                                     dic['rmin'][idx]),
                                 key = lambda t: t[1]))
    else:
      x_label = whichx
      idx, x_using = zip(*sorted(zip(idx,
                                     np.log10(dic[whichx][idx])),
                                 key = lambda t: t[1]))
      
    idx = list(idx)
    minx = np.min(x_using)
    maxx = np.max(x_using)
    
    if xRange == None:
        xRange = (minx, maxx)
        
    if figsize == None:
        figsize = (10, 8)
        
    matplotlib.rcParams['axes.linewidth'] = 1

    for item_group in lis:
        # All the items in a group in the list are plotted in the same pannel
        if ax == None:
            fig = plt.figure(figsize=figsize)
        
        if 'vr' in item_group[0]:
            yRange = item_group[0]['vr']
        else:
            val_tmp = dic[item_group[0]['name']][idx]
            vmax = np.max(val_tmp)
            vmin = np.min(val_tmp)
            yRange = (max(vmin/1.5, vmax/1e10), vmax * 4)
        
        if ax == None:
            ax = plt.subplot2grid((1, 1), (0, 0), colspan=1,
                              xlabel=x_label, ylabel='',
                              autoscalex_on=False, autoscaley_on=False,
                              xscale=xscale, yscale=yscale,
                              xlim=xRange, ylim=yRange)
        if not noTitle:
          plt.text(0.5, 0.98, 'z = {0:.2f} AU'.format(zpos), transform = plt.gca().transAxes,
             fontsize=25, horizontalalignment='center', verticalalignment='top')
        if set_ax_format:
            set_axis_format(ax, xscale=xscale, yscale=yscale)

        icolor = 0
        for item in item_group:
            name = item['name']
            if color == None:
                thiscolor = color_list[icolor%len(color_list)]
            else:
                thiscolor = color
            if label == None:
                if display_dic == None:
                  thislabel = name
                else:
                  thislabel = display_dic[name]
            else:
                thislabel = label
            if do_interpol:
                if dosmooth:
                    if winwidth == None:
                        winwidth=interpolNum/30
                else:
                    winwidth = 0
                xplot, yplot = loglin_interpol(x_using, dic[name][idx],
                                         xRange=xRange, num=interpolNum,
                                         method='linear', dosmooth=dosmooth, winwidth=winwidth)
                ax.plot(xplot[winwidth/2 : -1-winwidth/2],
                        yplot[winwidth/2 : -1-winwidth/2],
                        label=thislabel,
                        linewidth=lw, linestyle=linestyle,
                        color=thiscolor)
            else:
                ax.plot(x_using, dic[name][idx], label=thislabel,
                        linewidth=lw, linestyle=linestyle, marker='o', markersize=markersize,
                        color=thiscolor)
            icolor += 1
        
        if not noLegend:
          lgd = ax.legend(loc='best', prop={'size':15},
                          fancybox=False, shadow=False, ncol=1)
          lgd.get_frame().set_alpha(0.5)
    if returndata:
      return (x_using, dic[name][idx])
    elif returnaxis:
      return ax


def calc_colden(d, name, is_abundance=True):
    from numpy import zeros
    colden_name = 'N(' + name + ')'
    if is_abundance:
        N_cell = d['n_gas'] * d[name] * (d['zmax'] - d['zmin']) * 1.5E13
    else:
        N_cell = d[name] * (d['zmax'] - d['zmin']) * 1.5E13
    zmin = d['zmin']
    nc = len(zmin)
    N = zeros(nc)
    for i in range(nc):
        if i==0 or zmin[i] > zmin[i-1]:
            N[i] = N_cell[i]
        else:
            N[i] = N[i-1] + N_cell[i]
    d[colden_name] = N
    return



def make_fname_from_namelist(name_list):
    s = ''
    chars_to_avoid = '#%&{}\\<>*?/ &^$!`\'":@+|='
    for item in name_list:
        for it in item:
            s = s + '_' + it['name']
    t = ''
    for c in s[1:]:
        if c not in chars_to_avoid:
            t = t + c
        else:
            t = t + '_'
    return t


def overplot(ds, xr=(1,2), yr=(1e-12, 3e-4), zv=0.0, name='H2O', ylab='', linestyle='-',
             noTitle=True, ax=None, noLegend=False):
    import matplotlib.pyplot as plt
    for i,d in enumerate(ds):
        if ax == None:
            ax = plot_model_results_radial(d, [[{'name':name, 'vr':yr}]], zv, noTitle=noTitle, returnaxis=True,
                               color=mycolors[i], label='Iter'+str(i), do_interpol=True, interpolNum=500,
                               linestyle=linestyle,
                               xRange=xr, returndata=False, noLegend=True, figsize=(10,6))
        else:
            ax = plot_model_results_radial(d, [[{'name':name, 'vr':yr}]], zv, noTitle=noTitle, returnaxis=True,
                               color=mycolors[i], label='Iter'+str(i), ax=ax, do_interpol=True, interpolNum=500,
                               linestyle=linestyle,
                               xRange=xr, returndata=False, noLegend=True)
    ax.xaxis.grid(which='minor', color=(0.5,0.9,0.7), linewidth=0, linestyle='none')
    ax.yaxis.grid(which='minor', color=(0.5,0.9,0.7), linewidth=0, linestyle='none')
    ax.set_ylabel(ylab)
    if not noLegend:
        plt.legend(prop={'size':20}, framealpha=0.4, columnspacing=1, handletextpad=0)
    return ax



def stepplot(x1, x2, y, minval=0.0, color='blue', linewidth=5):
    import matplotlib.pyplot as plt
    for i in range(len(x1)):
        plt.plot([x1[i], x1[i], x2[i], x2[i]],
             [minval, y[i], y[i], minval], color=color, linewidth=linewidth)



def Bernstein(n, t):
    """For the generation of Bezier curve.
    n: degree
    t: 0 <= t <= 1
    The output is a 1-d array with n+1 elements.
    """
    import numpy as np
    
    if n == 0:
        return np.array([1])
    if n == 1:
        return np.array([1-t, t])
    
    b = np.array([1-t, t])
    
    for j in range(1,n):
        btmp_1 = b*(1-t)
        btmp_2 = b*t
        b = np.concatenate((btmp_1[0:1], btmp_1[1:(j+1)]+btmp_2[0:(j)], btmp_2[j:]))
    
    return b


def Bezier_Curve(xin, yin, n_points):
    """
    - Input
      1. X: All control points.
            Its format must be:
              [[x1, x2, ..., xn], $
               [y1, y2, ..., yn]]
    
      2. n_points: Number of points in the output curve.
    
    - Output
      1. Xn: Coordinates of the output curve,
             in the same format of X.
    
     1. No validation check are made.
     2. It is easy to extend this algorithm to higher dimensions.
    """
    import numpy as np
    
    n = xin.size - 1
    
    t = np.linspace(0, 1, n_points)
    
    xout = np.zeros(n_points)
    yout = np.zeros(n_points)
    
    for i in range(0, n_points):
        btmp = Bernstein(n, t[i])
        xout[i] = np.sum(xin*btmp)
        yout[i] = np.sum(yin*btmp)
    
    return (xout, yout)



def plot_a_pannel_bezier(isub, dic, idx, name, xin, xr, yr, minval=1e-50,
                         xscale='linear', yscale='log', npt_bezier=300):
    import matplotlib.pyplot as plt
    plt.subplot(*isub); isub[2] += 1
    yin = dic[name][idx]
    xout, yout = Bezier_Curve(xin, yin, npt_bezier)

    plot(xout, yout, linewidth=2)
    gca().set_xlim(xr)
    gca().set_ylim(yr)
    gca().set_xscale(xscale)
    gca().set_yscale(yscale)
    set_axis_format(gca(), xscale=xscale, yscale=yscale)
    text(0.5, 0.97, name, transform = gca().transAxes, fontsize=15,
         horizontalalignment='center', verticalalignment='top', weight='bold')



colors_ = ['red', 'green', 'blue', 'CadetBlue', 'magenta', 'BurlyWood',
           'Brown', 'DarkOliveGreen', 'DarkSlateBlue', 'DodgerBlue', 'BlueViolet', 'Chocolate',
           'Tomato', 'Olive', 'SteelBlue', 'SlateBlue', 'Orchid', 'Orange', 'black']
linestyles_ = ['-', '--', ':']
linestyles = []
for lstyle_ in linestyles_:
    linestyles += [lstyle_] * len(colors_)
mycolors = colors_ * len(linestyles_)
linewidths = [4] * len(mycolors)
del colors_, linestyles_


def plot_a_pannel_simple(x, y, name, isub=[1,1,1], xr=None, yr=None,
                  xscale='log', yscale='log'):
    import numpy as np
    import matplotlib.pyplot as plt
    plt.subplot(*isub); isub[2] += 1
    if xr == None:
        xr = (np.min(x), np.max(x))
    if yr == None:
        idx = np.where(np.logical_and(x >= xr[0], x <= xr[1]))
        yr = (np.min(y[idx])/2, 5*np.max(y[idx]))
    plt.plot(x, y, linewidth=3)
    plt.gca().set_xlim(xr)
    plt.gca().set_ylim(yr)
    plt.gca().set_xscale(xscale)
    plt.gca().set_yscale(yscale)
    set_axis_format(plt.gca(), xscale=xscale, yscale=yscale)
    plt.text(0.5, 0.97, name, transform = plt.gca().transAxes, fontsize=15, weight='bold',
         horizontalalignment='center', verticalalignment='top')
    
    
def plot_a_pannel_simple_many(xylist, name, isub=None, xr=None, yr=None,
                  xscale='log', yscale='log', majorgridon=True, minorgridon=True):
    colors = ['red', 'green', 'blue', 'CadetBlue', 'magenta', 'BurlyWood',
              'Brown', 'DarkOliveGreen', 'DarkSlateBlue', 'DodgerBlue', 'BlueViolet', 'Chocolate',
              'Tomato', 'Olive', 'SteelBlue', 'SlateBlue', 'Orchid', 'Orange', 'black']

    import numpy as np
    import matplotlib.pyplot as plt
    if isub == None:
      isub = [1,1,1]
    plt.subplot(*isub); isub[2] += 1
    
    for i in range(len(xylist)):
        lw = 2*(3-2*(i%2))
        alph = (i%2)*0.5 + 0.5
        if len(xylist[i]) > 2:
            plt.plot(xylist[i][0], xylist[i][1], alpha=alph,
                 linewidth=lw, label=xylist[i][2], color=colors[i])
        else:
            plt.plot(xylist[i][0], xylist[i][1], alpha=alph,
                 linewidth=lw, color=colors[i])
    if len(xylist[0]) > 2:
        lgd = plt.gca().legend(loc='best', prop={'size':15},
                           fancybox=False, shadow=False, ncol=1)
        lgd.get_frame().set_alpha(0.2)

    plt.gca().set_xlim(xr)
    plt.gca().set_ylim(yr)
    plt.gca().set_xscale(xscale)
    plt.gca().set_yscale(yscale)
    set_axis_format(plt.gca(), xscale=xscale, yscale=yscale)
    plt.text(0.5, 0.97, name, transform = plt.gca().transAxes, fontsize=15, weight='bold',
         horizontalalignment='center', verticalalignment='top')


def quiver_plot(d, zname, uname, vname, 
                fsize = (16,10),
                xr = (0, 4), yr = (0, 2),
                vmin=0.0, vmax=1.0,
                xscale='linear', yscale='linear',
                xlabel='r (AU)', ylabel='z (AU)'):
    import matplotlib.pyplot as plt
    X = 0.5 * (d['rmin'] + d['rmax'])
    Y = 0.5 * (d['zmin'] + d['zmax'])
    Z = d[zname]
    U = d[uname]
    V = d[vname]
    xi = linspace(xr[0], xr[1], num=100)
    yi = linspace(yr[0], yr[1], num=100)
    ZZ = griddata(X, Y, Z, xi, yi)
    figure(figsize=fsize)
    ax = plt.subplot2grid((1, 15), (0, 0), colspan=14,
                          xlabel=xlabel, ylabel=ylabel,
                          autoscalex_on=False, autoscaley_on=False,
                          xscale=xscale, yscale=yscale,
                          xlim=xr, ylim=yr)
    C = plt.contourf(xi, yi, ZZ,
                     cmap=plt.cm.rainbow,
                     vmin=vmin, vmax=vmax,
                     extend='neither')
    #C.cmap.set_under('black')
    #C.cmap.set_under('white')
    plt.colorbar(C)
    ax.quiver(X, Y, U, V)
    set_axis_format(ax)

def colorize_axis(fig, _ax, color='blue', alpha=0.5, pad=0.0):
    from matplotlib.transforms import Bbox
    from matplotlib.patches import Rectangle
    _ax.figure.canvas.draw()
    items = _ax.get_xticklabels() + _ax.get_yticklabels()
    items += [_ax, _ax.title]
    bbox = Bbox.union([item.get_window_extent() for item in items])
    extent = bbox.expanded(1.0 + pad, 1.0 + pad)
    extent = extent.transformed(fig.transFigure.inverted())
    rect = Rectangle([extent.xmin, extent.ymin], extent.width, extent.height,
                     facecolor=color, edgecolor='none', alpha=alpha,
                     zorder=_ax.get_zorder()-1, 
                     transform=fig.transFigure)
    fig.patches.append(rect)

def get_val(d, r, z, k):
    for i in range(len(d['rmin'])):
        if d['rmin'][i] <= r and r <= d['rmax'][i] \
            and d['zmin'][i] <= z and z <= d['zmax'][i]:
            return k, r, z, d[k][i]


def rgb2brightness(r, g, b):
    # http://www.nbdtech.com/Blog/archive/2008/04/27/Calculating-the-Perceived-Brightness-of-a-Color.aspx
    import numpy as np
    return np.sqrt(0.241*r*r + 0.691*g*g + 0.068*b*b)

def plot_cm_brightness(cm):
    import numpy as np
    import matplotlib.pyplot as plt
    x = np.arange(0, 1, 0.005)
    _c = cm(x)
    y = np.array([rgb2brightness(*_[0:3]) for _ in _c])
    plt.figure(figsize=(10,10))
    plt.scatter(x, y, c=x, cmap=cm, edgecolor='none', s=100, marker=',')
    plt.plot(x, y, color='white')
    set_axis_format(plt.gca())
    plt.tight_layout()

def make_my_colormap(c_list=None):
    if c_list == None:
        c_list = [(0.0, 'blue'),
                  (0.2, 'CornflowerBlue'),
                  (0.4, 'cyan'),
                  (0.6, 'Chartreuse'),
                  (0.75, 'yellow'),
                  (1.0, 'red')]
    import matplotlib.colors as colors
    return colors.LinearSegmentedColormap.from_list('mycm', c_list)


def zoomin_annot(fig, ax, rect, xy1, xy2, edgecolor='Brown', linewidth=1, linecolor='Brown', zorder=100, clip_on=False):
    import matplotlib as mpl
    from matplotlib.path import Path
    _x = rect[0]; _y = rect[1]; _dx = rect[2]; _dy = rect[3]
    ax.add_patch(mpl.patches.Rectangle((_x,_y), _dx, _dy, edgecolor=edgecolor,
                                       facecolor='none',
                                       zorder=100, transform=ax.transAxes))
    _verts = [fig.transFigure.inverted().transform(tuple(ax.transAxes.transform((_x, _y+_dy)))),
              xy1,
              fig.transFigure.inverted().transform(tuple(ax.transAxes.transform((_x+_dx,_y+_dy)))),
              xy2]
    _codes = [Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO]
    _path = Path(_verts, _codes)
    ax.add_patch(mpl.patches.PathPatch(_path, lw=linewidth, color=linecolor, clip_on=clip_on,
                                       transform=fig.transFigure, zorder=zorder))


def get_v_thermal(T, m):
    import numpy as np
    return np.sqrt(8.0*1.38e-16*T / (np.pi * m * 1.67e-24))

def get_Tevap(E_ev, n_H, Tgas=1e2, d2gmass=1e-2, rd=0.1, rhod=2.0, molmass=20.0, nu=1e12, freeze_threshold=1e1):
    import numpy as np
    v_T = get_v_thermal(Tgas, molmass)
    return E_ev / (np.log(freeze_threshold / (v_T/nu * n_H * (3.0*1.67e-24*d2gmass/(4.0*(rd*1e-4)*rhod)))
                   ))



def colorful_hist(x, y, w, sepy, bins=20, iflog=False, xr=None, cm_=None):
    import numpy as np
    import matplotlib.pyplot as plt
    if cm_ == None:
        import matplotlib.cm as cm
        cmap = cm.rainbow
    else:
        cmap = cm_
    idx = np.digitize(y, sepy)
    idxu = np.sort(np.unique(idx))
    nu = idxu.size
    nsep = len(sepy)
    indices = [idx==idxu[i] for i in range(nu)]
    hlist = [x[indices[i]] for i in range(nu)]
    wlist = [w[indices[i]] for i in range(nu)]
    clist = [cmap(float(idxu[i])/(nsep+1)) for i in range(nu)]
    h = plt.hist(hlist, bins=bins, weights=wlist, color=clist,
             log=iflog, range=xr, edgecolor='none', linewidth=0,
             histtype='barstacked', rwidth=0.95,
             orientation='vertical', stacked=True, antialiased=True, fill=True)



def data2axes(fig, ax, xy):
    """Convert data coordinate into axes coordinate."""
    import numpy as np
    return fig.transFigure.inverted().transform(ax.transData.transform(np.array([[xy[0], xy[1]]])))

def data2axes1(ax, xy):
    """Convert data coordinate into axes coordinate."""
    import numpy as np
    return ax.transAxes.inverted().transform(ax.transData.transform(np.array([[xy[0], xy[1]]])))



def get_col_integ(dic, name_x, name_y, name_sel, sel_range,
                  f_filter=None, param_filter=None,
                  do_reg=False, fill_blank=False,
                  n=50, dx=None, xr=None, dx_tol=1e-5, winw=2.0, y_thresh=1e-1):
    import numpy as np
    AU2cm = 1.5e13
    all_x = dic[name_x]
    all_y = dic[name_y]
    all_s = dic[name_sel]
    #
    uni_x, iuni = np.unique(all_x.round(decimals=5), return_inverse=True)
    iuni = np.array(iuni)
    #
    uni_y = np.zeros_like(uni_x)
    for i,x in enumerate(uni_x):
        idx = np.where(np.logical_and(i == iuni,
                             np.logical_and(all_s >= sel_range[0],
                                            all_s <= sel_range[1])))[0]
        if f_filter != None and param_filter != None:
            idx = idx[f_filter(idx, **param_filter)]
        uni_y[i] = np.sum((dic['zmax'][idx]-dic['zmin'][idx])
                          * dic['n_gas'][idx]
                          * all_y[idx]) * AU2cm * 2.0
    if fill_blank:
        for i in range(len(uni_x)):
            if uni_y[i] <= y_thresh:
                if i > 0 and i < len(uni_x)-1:
                    fL = False
                    fR = False
                    for j in range(i-1, -1, -1):
                        if uni_y[j] > y_thresh:
                            iL = j
                            fL = True
                            break
                    for j in range(i+1, len(uni_x), 1):
                        if uni_y[j] > y_thresh:
                            iR = j
                            fR = True
                            break
                    if fL and fR:
                        uni_y[i] = uni_y[iL] + \
                          (uni_y[iR] - uni_y[iL]) / (uni_x[iR] - uni_x[iL]) \
                          * (uni_x[i] - uni_x[iL])
    if not do_reg:
        return (uni_x, uni_y)
    else:
        minx = uni_x.min()
        maxx = uni_x.max()
        if xr != None:
            xmin = xr[0]
            xmax = xr[1]
        else:
            xmin = minx
            xmax = maxx
        if dx != None:
            x_reg = np.arange(xmin, xmax+dx, dx)
        else:
            x_reg = np.linspace(xmin, xmax, num=n)
        y_reg = np.zeros_like(x_reg)
        for i,x in enumerate(x_reg):
            dist = np.abs(x-uni_x)
            if x <= minx:
                winsize = uni_x[1] - uni_x[0]
            elif x >= maxx:
                winsize = uni_x[-1] - uni_x[-2]
            else:
                idx1 = np.logical_and(uni_x > x, uni_y > y_thresh)
                idx2 = np.logical_and(uni_x < x, uni_y > y_thresh)
                if np.any(idx1) and np.any(idx2):
                    winsize = np.min(uni_x[idx1]) - np.max(uni_x[idx2])
                else:
                    y_reg[i] = 0.0
                    continue
            winsize = (winsize + dx_tol) * winw
            idx = np.logical_and(x-winsize <= uni_x, x+winsize >= uni_x)
            if np.any(idx):
                y_reg[i] = np.average(uni_y[idx], weights=np.exp(-dist[idx]/(winsize/winw)))
            else:
                y_reg[i] = 0.0
        return (x_reg, y_reg)



def loglin_interpol(x, y, xRange=None, yRange=None, num=50, method='linear', dosmooth=False, winwidth=None, bounds_error=False):
    import numpy as np
    from scipy import interpolate
    #
    def whetherUseLinear(maxx, minx, threshold=1e2):
        if maxx*minx <= 0.0 or maxx < 0.0:
            uselinearX = True
        elif maxx/minx > threshold:
            uselinearX = False
        else:
            uselinearX = True
        return uselinearX
    #
    if xRange == None:
        maxx, minx = np.max(x), np.min(x)
    else:
        maxx, minx = xRange[0], xRange[1]
    if yRange == None:
        maxy, miny = np.max(y), np.min(y)
    else:
        maxy, miny = yRange[0], yRange[1]
    uselinearX = whetherUseLinear(maxx, minx)
    uselinearY = whetherUseLinear(maxy, miny)
    if uselinearX:
        xuse = x
        xnew = np.linspace(minx, maxx, num=num)
    else:
        xuse = np.log(x)
        xnew = np.linspace(np.log(minx), np.log(maxx), num=num)
    if uselinearY:
        yuse = y
    else:
        yuse = np.log(y)
    fint = interpolate.interp1d(xuse, yuse, kind=method, bounds_error=bounds_error, fill_value=np.nan)
    ynew = fint(xnew)
    if dosmooth:
      if winwidth == None or winwidth <= 2:
        winwidth = 2 + num/100
      window = np.ones(winwidth, 'd')
      ynew = np.convolve(ynew, window/window.sum(), mode='same')
    if not uselinearX:
        xnew = np.exp(xnew)
    if not uselinearY:
        ynew = np.exp(ynew)
    return xnew, ynew



def customaxis(ax, c_left='k', c_bottom='k', c_right='none', c_top='none',
               lw=1, size=12, pad=5, tickdir='in', ticklen=4):

    for c_spine, spine in zip([c_left, c_bottom, c_right, c_top],
                              ['left', 'bottom', 'right', 'top']):
        if c_spine != 'none':
            ax.spines[spine].set_color(c_spine)
            ax.spines[spine].set_linewidth(lw)
        else:
            ax.spines[spine].set_color('none')
    if (c_bottom == 'none') & (c_top == 'none'): # no bottom and no top
        ax.xaxis.set_ticks_position('none')
    elif (c_bottom != 'none') & (c_top != 'none'): # bottom and top
        ax.tick_params(axis='x', direction=tickdir, width=lw, length=ticklen,
                      color=c_bottom, labelsize=size, pad=pad)
    elif (c_bottom != 'none') & (c_top == 'none'): # bottom but not top
        ax.xaxis.set_ticks_position('bottom')
        ax.tick_params(axis='x', direction=tickdir, width=lw, length=ticklen,
                       color=c_bottom, labelsize=size, pad=pad)
    elif (c_bottom == 'none') & (c_top != 'none'): # no bottom but top
        ax.xaxis.set_ticks_position('top')
        ax.tick_params(axis='x', direction=tickdir, width=lw, length=ticklen,
                       color=c_top, labelsize=size, pad=pad)
    if (c_left == 'none') & (c_right == 'none'): # no left and no right
        ax.yaxis.set_ticks_position('none')
    elif (c_left != 'none') & (c_right != 'none'): # left and right
        ax.tick_params(axis='y', direction=tickdir, width=lw, length=ticklen,
                       color=c_left, labelsize=size, pad=pad)
    elif (c_left != 'none') & (c_right == 'none'): # left but not right
        ax.yaxis.set_ticks_position('left')
        ax.tick_params(axis='y', direction=tickdir, width=lw, length=ticklen,
                       color=c_left, labelsize=size, pad=pad)
    elif (c_left == 'none') & (c_right != 'none'): # no left but right
        ax.yaxis.set_ticks_position('right')
        ax.tick_params(axis='y', direction=tickdir, width=lw, length=ticklen,
                       color=c_right, labelsize=size, pad=pad)


def get_warm_water_mass(c, Twarm=1.5e2):
    name = 'H2O'
    idx = where(c['Tgas'] > Twarm)
    mass = sum(c['mg_cell'][idx] * c[name][idx] * 18.0 / 1.4) * 2
    Earth_ocean_mass = 1.4e24 # gram
    return mass / Earth_ocean_mass

def get_warm_water_volume(c):
    name = 'H2O'
    idx = where(c['Tgas'] > 1.5e2)
    def get_volume(rmin, rmax, zmin, zmax):
        return pi * (rmax**2 - rmin**2) * (zmax - zmin) * 2
    v = 0.0
    for i in idx[0]:
        v += get_volume(c['rmin'][i], c['rmax'][i], c['zmin'][i], c['zmax'][i])
    return v

def get_diffuse_water_mass(c, Tcold=1e2):
    name = 'H2O'
    idx = where(c['Tdust'] < Tcold)
    mass = sum(c['mg_cell'][idx] * c[name][idx] * 18.0 / 1.4) * 2
    Earth_ocean_mass = 1.4e24 # gram
    return mass / Earth_ocean_mass

def get_all_water_mass(c):
    name = 'H2O'
    mass = sum(c['mg_cell'] * c[name] * 18.0 / 1.4) * 2
    Earth_ocean_mass = 1.4e24 # gram
    return mass / Earth_ocean_mass

def get_medianG0_diffuse_water_mass(c):
    name = 'H2O'
    idx = where(c['Tdust'] < 1e2)
    return median(c['G0_UV'][idx])

def plot_a_vs_b(x_, y_, dic=None, xr=None, yr=None,
                c=None, clog=False, crange=None, cmap=None,
                s=None, slog=False, srange=None,
                xs='log', ys='log',
                linthreshx=1e-5, linthreshy=1e-5,
                f_=None, ax_=None, pos_=None,
                figsize=(10,6)):
    import numpy as np
    import matplotlib.pyplot as plt
    if dic != None:
        x = dic[x_]
        y = dic[y_]
    else:
        x = x_
        y = y_
    if c != None:
        if crange != None:
            c_ = np.clip(c, crange[0], crange[1])
        else:
            c_ = c
        if clog:
            c_ = np.log10(c_)
    else:
        c_ = 'b'
    if s != None:
        if srange != None:
            s_ = np.clip(s, srange[0], srange[1])
        else:
            s_ = s
        if slog:
            s_ = np.log(s_)
        smax = s_.max()
        smin = s_.min()
        min_s = 10
        range_s = 200
        s_ = min_s + (s_ - smin) * (range_s / (smax - smin))
    else:
        s_ = 20
    if f_ == None:
        f = plt.figure(figsize=figsize)
    else:
        f = f_
    if ax_ == None:
        if pos_ == None:
            pos = (0.12,0.15,0.7,0.8)
        else:
            pos = pos_
        ax = f.add_axes(pos, xlim=xr, ylim=yr,
                        xscale=xs, yscale=ys)
        ax.set_xscale(xs, linthreshx=linthreshx)
        ax.set_yscale(ys, linthreshy=linthreshy)
    else:
        ax = ax_
    C = ax.scatter(x, y, edgecolor='none', c=c_, cmap=cmap, s=s_, alpha=0.7)
    if c != None:
        plt.colorbar(C)
    set_axis_format(ax, xscale=xs, yscale=ys)
    return f, ax




hlist = [\
    {'h_ph_gr': 'Photoelec'},
    {'h_fo_H2': 'H$_2$ form'},
    {'h_cosmi': 'Cosmic-ray'},
    {'h_vi_H2': 'H$_2$ vib'},
    {'h_io_CI': r'C $\textsc{i}$'},
    {'h_ph_H2': 'H$_2$ photodiss'},
    {'h_ph_wa': 'Water photodiss'},
    {'h_ph_OH': 'OH photodiss'},
    {'h_Xray' : 'X-ray'},
    {'h_visco': 'Viscosity'},
    {'h_chem' : 'Chemical'},
    {'h_gg_co': 'Gas-grain colli'},]

clist =[\
    {'c_el_gr': 'Photoelec'},
    {'c_vi_H2': 'H$_2$ vib'},
    {'c_chem' : 'Chemical'},
    {'c_gg_co': 'Gas-grain colli'},
    {'c_OI'   : r'O $\textsc{i}$'},
    {'c_CII'  : r'C $\textsc{ii}$'},
    {'c_NII'  : r'N $\textsc{ii}$'},
    {'c_SiII'  : r'Si $\textsc{ii}$'},
    {'c_FeII'  : r'Fe $\textsc{ii}$'},
    {'c_OH_ro': 'OH rotat'},
    {'c_wa_ro': 'Water rotat'},
    {'c_wa_vi': 'Water vib'},
    {'c_CO_ro': 'CO rotat'},
    {'c_CO_vi': 'CO vib'},
    {'c_H2_ro': 'H$_2$ rotat'},
    {'c_LyAlp': r'Lyman $\alpha$'},
    {'c_fb'   : 'Free-free'},
    {'c_ff'   : 'Free-bound'},]

display_dic = {'H2O': 'H$_2$O',
               'gH2O': 'H$_2$O ice',
               'H2': 'H$_2$',
               'O2': 'O$_2$',
               'N_H2O_I': '$N_{\\rm water, toISM}$',
               'N_CO_I': '$N_{\\rm CO, toISM}$',
               'N_H2O_S': '$N_{\\rm water, toStar}$',
               'N_CO_S': '$N_{\\rm CO, toStar}$',
               'N_H2_S': '$N_{\\rm H2, toStar}$',
               'Tgas': '$T_{\\rm gas}$',
               'n_gas':'$n_{\\rm gas}$',
               'Tdust': '$T_{\\rm dust}$',
               'ndust_t':'$n_{\\rm dust}$',
               'LyAG0_a': '$G_{0,{\\rm Ly}\\alpha}$',
               'G0_UV': '$G_{0,{\\rm UV}}$',
               'flx_NIR': 'NIR flux',
               'presr_g': r'$P_g$',
               'presr_t': r'$P_t$',
               'pres_t2g': r'$P_{t2g}$',
               'h_visco': r'$h_{\rm visco}$',
               'Ncol_I': r'$N_{\rm ISM}$',
               'Av_Star': r'Av$_{\rm Star}$',
               'Av_ISM': r'Av$_{\rm ISM}$',
               'flx_Lya': r'flux$_{\rm Lya}$',
               'flx_Vis': r'flux$_{\rm vis}$',
               'flx_Xray': r'flux$_{\rm Xray}$',
               'flx_NIR': r'flux$_{\rm NIR}$',
               'f_H2_I': r'$f_{\rm shield, I}$ (H$_2$)',
               'f_H2_S': r'$f_{\rm shield, S}$ (H$_2$)',
               'f_CO_I': r'$f_{\rm shield, I}$ (CO)',
               'f_CO_S': r'$f_{\rm shield, S}$ (CO)',
               'zeta_X': r'$\zeta_{\rm X}$',
               'h_Xray' : r'h$_{\rm X}$',
               }


def smooth_spec(f, spec, f0, df0, resolvingpower):
    from numpy import floor, zeros, convolve, zeros_like, arange, exp, log
    c = 2.99792458e8
    df_spectrometer = f0 / resolvingpower
    specSmooth = floor(df_spectrometer / df0)
    #
    n = len(f)
    specExtend = specSmooth*4
    stmp = zeros(specExtend*2+n)
    stmp[0:specExtend] = spec[0]
    stmp[specExtend:specExtend+n] = spec
    stmp[specExtend+n:] = spec[-1]
    ftmp = zeros_like(stmp)
    ftmp[specExtend:specExtend+n] = f
    ftmp[0:specExtend] = f[0] - arange(specExtend, 0, -1) * (f[1] - f[0])
    ftmp[specExtend+n:] = f[-1] + arange(1, specExtend+1, 1) * (f[1] - f[0])
    #
    cv = arange(specSmooth*4) - specSmooth*2
    cv = exp(-log(16.0)/specSmooth**2 * cv**2)
    #cv = exp(-(cv/(specSmooth*0.5))**2)
    #cv = zeros(specSmooth*4)
    #cv[specSmooth : specSmooth*2] = 1.0
    cv /= sum(cv)
    sp = convolve(stmp, cv, mode='same')
    lam = c/ftmp * 1e6
    sp = sp[specSmooth*1.5:-specSmooth*1.5]
    lam  = lam[specSmooth*1.5:-specSmooth*1.5]
    return lam, sp


def plot_spectrum_combine(h, lw=1, visible=True, printout=True,
                  axspec=None, figspec=None, figspecsize=(7,4),
                  lamRange=None, ySpecRange=None, lgd_str='',
                  xscale='linear', yscale='linear', unit='Jy',
                  removebase=False, showlegend=False, color='black',
                  xtype=None, vel_range=None, beamsize=None,
                  resolvingpower=None, qnum_select=None,
                  returndata=False, set_ax=True,
                  returnax=False):
    from numpy import arange
    n_ext = len(h)
    c = 2.99792458e8
    for ext in h:
        hd = ext.header
        #print hd
        if hd['EXTNAME'] == 'FluxSpec':
            spec = ext.data[0]
            n = hd['NAXIS1']
            k0 = hd['CRPIX1']
            f0 = hd['F0']
            E_up = hd['EUP']
            if 'QNUM' in hd:
                qnum = hd['QNUM']
            else:
                qnum = ''
            Bul = hd['BUL']
            if hd['CTYPE1'] == 'V':
                dv = hd['CDELT1']
                v0 = hd['CRVAL1']
                v = (arange(1, n+1) - k0) * dv + v0
                f = f0 * (1-v/c)
                df = abs(f[1]-f[0])
            elif hd['CTYPE1'] == 'F':
                fmin = hd['CRVAL1']
                df = hd['CDELT1']
                f = (arange(1, n+1) - k0) * df + fmin
            lam0 = hd['LAM0']/1e4
            lam = c/f*1e6
            break

    if resolvingpower != None:
        #df_spectrometer = f[len(f)/2] / resolvingpower
        #specSmooth = floor(df_spectrometer / abs(f[1]-f[0]))
        ##
        #specExtend = specSmooth*4
        #stmp = zeros(specExtend*2+n)
        #stmp[0:specExtend] = spec[0]
        #stmp[specExtend:specExtend+n] = spec
        #stmp[specExtend+n:] = spec[-1]
        #ftmp = zeros_like(stmp)
        #ftmp[specExtend:specExtend+n] = f
        #ftmp[0:specExtend] = f[0] - arange(specExtend, 0, -1) * (f[1] - f[0])
        #ftmp[specExtend+n:] = f[-1] + arange(1, specExtend+1, 1) * (f[1] - f[0])
        ##
        ##cv = array([1.0]*specSmooth)
        #cv = arange(specSmooth*4) - specSmooth*2
        #cv = exp(-(cv/(specSmooth*0.5))**2)
        #cv /= sum(cv)
        #spec = convolve(stmp, cv, mode='same')
        #lam = c/ftmp * 1e6
        #spec = spec[specSmooth*1.5:-specSmooth*1.5]
        #lam  = lam[specSmooth*1.5:-specSmooth*1.5]
        lam, spec = smooth_spec(f, spec, f0, abs(f[1]-f[0]), resolvingpower)
    if lamRange == None:
        lamRange = (lam.min(), lam.max())
    if ySpecRange == None:
        yrange_ = spec.max() - spec.min()
        ySpecRange = (spec.min() - yrange_*0.1, spec.max() + yrange_*0.15)

    if figspec == None:
        figspec = figure(figsize=figspecsize)
    if axspec == None:
        if unit == 'Jy':
            ylabel = 'Flux (Jy)'
        if unit == 'mJy':
            ylabel = 'Flux (mJy)'
        if unit == 'nuFnu':
            ylabel = r'$\nu F_{\nu}$ (erg s$^{-1}$ cm$^{-2}$)'
        if xtype == 'vel':
            axspec = figspec.add_axes((0.1,0.1,0.85,0.85),
                              xlabel=r'v (km s$^{-1}$)',
                              ylabel=ylabel,
                              autoscalex_on=False, autoscaley_on=False,
                              xscale=xscale, yscale=yscale,
                              xlim=vel_range, ylim=ySpecRange)
        else:
            axspec = figspec.add_axes((0.1,0.1,0.85,0.85),
                              xlabel=r'$\lambda$ ($\mu$m)',
                              ylabel=ylabel,
                              autoscalex_on=False, autoscaley_on=False,
                              xscale=xscale, yscale=yscale,
                              xlim=lamRange, ylim=ySpecRange)
    if unit == 'mJy':
        spec *= 1e3
    if unit == 'nuFnu':
        f = 2.99792458e8 / (lam*1e-6)
        spec = spec * (f*1e-23)
    if removebase:
        a = (spec[-1] - spec[0]) / (lam[-1] - lam[0])
        b = spec[0]
        spec = spec - (a * lam - (a * lam[0] - b))
    #
    fmtstr = '{:s}={:8.3f}, {:s}={:12.4e}, {:s}={:6.0f}, {:s}={:s}, {:s}={:9.2e}, {:s}={:10.4f}'
    if qnum_select != None:
      if not (qnum_select in qnum):
        fmtstr = ''
    if xtype == 'vel':
        vel = (lam/lam0-1.0)*c/1e3
        axspec.plot(vel, spec, linewidth=lw, label=lgd_str, visible=visible, color=color)
        #
        if visible != False and printout and fmtstr != '':
          if unit == 'Jy' and removebase:
              print(fmtstr.format(\
                              'Lam', lam0,
                              'Freq', f0,
                              'Eup', E_up,
                              'Qnum', qnum, 
                              'Bul', Bul, 
                              'IntFdnu (1e-18 W m-2)', sum(spec) * 1e-26 * df * 1e18))
              if beamsize != None:
                k = 1.38e-23
                Omega = theta2Omega( beamsize * pi / (3600*180.0) )
                print('IntTdV (mK km/s): ', sum(spec) * 1e-26 * df * (lam0*1e-6)**3 / (2*k*Omega) / 1e3 * 1e3)
          if unit == 'mJy' and removebase:
              print(fmtstr.format(\
                              'Lam', lam0,
                              'Freq', f0,
                              'Eup', E_up,
                              'Qnum', qnum, 
                              'Bul', Bul, 
                              'IntFdnu (1e-18 W m-2)', sum(spec) * 1e-26 * df * 1e18 * 1e-3))
              if beamsize != None:
                k = 1.38e-23
                Omega = theta2Omega( beamsize * pi / (3600*180.0) )
                print('IntTdV (mK km/s): ', sum(spec) * 1e-26 * df * (lam0*1e-6)**3 / (2*k*Omega) / 1e3)
    else:
        axspec.plot(lam, spec, linewidth=lw, label=lgd_str, visible=visible, color=color)
        #
        if visible != False and printout and fmtstr != '':
          if unit == 'Jy':
              print(fmtstr.format(\
                              'Lam', lam0,
                              'Freq', f0,
                              'Eup', E_up,
                              'Qnum', qnum, 
                              'Bul', Bul, 
                              'IntFdnu (1e-18 W m-2)', sum(spec) * 1e-26 * df * 1e18))
          if unit == 'mJy':
              print(fmtstr.format(\
                              'Lam', lam0,
                              'Freq', f0,
                              'Eup', E_up,
                              'Qnum', qnum, 
                              'Bul', Bul, 
                              'IntFdnu (1e-18 W m-2)', sum(spec) * 1e-26 * df * 1e18 * 1e-3))
    if set_ax:
        set_axis_format(axspec, color_major=(0.5,0.5,0.5), color_minor=(0.8,0.8,0.8))
    if showlegend:
        axspec.legend(fontsize=20, framealpha=0.6, columnspacing=1, handletextpad=0)
    if returnax and not returndata:
        return axspec, figspec
    if returndata and not returnax:
        return lam, spec
    if returnax and returndata:
        return axspec, figspec, lam, spec


def TpeakFWHM2IntT(T, FWHM):
    '''
    T: K
    FWHM: km/s
    Return the integrated intensity in K km/s
    '''
    return T * FWHM / sqrt(4*log(2)/pi)

def IntT2Tpeak(IntT, FWHM):
    '''
    IntT: in K km/s
    FWHM: in km/s
    Return the peak temperature in K.
    '''
    return IntT/FWHM * sqrt(4*log(2)/pi)

def T2Flux(T, nu, Omega):
    '''
    T: in K
    nu: in Hz
    Omega: beam size in steradian
    Return the flux in SI unit.
    '''
    c = 2.99792458e8
    k = 1.38e-23
    lam = c/nu
    return 2*k*T/lam**2 * Omega

def IntT2Flux(IntT, FWHM, nu, Omega):
    '''
    IntT: in K km/s
    FWHM: in km/s
    nu: in Hz
    Omega: beam size in steradian
    Return the peak flux in SI unit.
    '''
    return T2Flux(IntT2Tpeak(IntT, FWHM), nu, Omega)

def IntT2IntFlux(IntT, nu, Omega):
    '''
    IntT: in K km/s
    nu: in Hz
    Omega: beam size in steradian
    Return the integrated flux in SI unit.
    '''
    k = 1.38e-23
    c = 2.99792458e8
    lam = c/nu
    return 2*k*Omega/lam**3 * (IntT * 1e3)

def theta2Omega(theta):
    '''
    theta is the FWHM of the main beam.
    Return the integrated beam size.
    '''
    return pi / log(16) * theta**2


def jy2nuFnu(spec, lam):
    ''' spec in Jy
        lam in micron
    '''
    f = 2.99792458e8 / (lam*1e-6)
    return spec * (f*1e-23)




def coadd_spec(d, dx=None, xmin=None, xmax=None):
    import scipy.interpolate
    from numpy import linspace, zeros_like, interp
    if xmin == None:
        xmin = min([min(_[0][0], _[0][-1]) for _ in d])
    if xmax == None:
        xmax = max([max(_[0][0], _[0][-1]) for _ in d])
    if dx == None:
        #dx = min([abs(_[0][1]-_[0][0]) for _ in d])
        dx = (xmax - xmin) / 1e4
    x = linspace(xmin, xmax, num=int((xmax-xmin)/dx)+1)
    y = zeros_like(x)
    dx = x[1]-x[0]
    for t in d:
        if t[0][0] > t[0][-1]:
            xx = t[0][-1::-1]
            yy = t[1][-1::-1]
        else:
            xx = t[0]
            yy = t[1]
        i0 = int((xx[0]-xmin)/dx)
        i1 = int((xx[-1]-xmin)/dx)+1
        y[i0:i1] += interp(x[i0:i1], xx, yy, left=0.0, right=0.0)
    return x, y


def plot_spectra_dir(files, lamRange=(27, 35), yRange=(-0.1, 0.6), ax=None, fig=None,
                     color='black', color_coadd='red', hide_comp=False,
                     lw=1.0, lw_coadd=1.0,
                     returndata=False, set_ax=True,
                     resolvingpower=6e2, coadd=True, removebase=True, printout=True, qnum_select=None):
    from astropy.io import fits
    import matplotlib.pyplot as plt
    lamRange = lamRange
    resolvingpower = resolvingpower
    if ax == None or fig == None:
        filename = files[0]
        h3 = fits.open(filename)
        ax, fig = plot_spectrum_combine(h3, lamRange=lamRange, ySpecRange=yRange, unit='Jy', set_ax=set_ax,
                                        xscale='linear', yscale='linear', visible=False, printout=printout,
                                        figspecsize=(16,4), returnax=True, lw=1)

    d = []
    for filename in files:
        h3 = fits.open(filename)
        ax, fig, lam, spec = plot_spectrum_combine(h3, returnax=True, returndata=True, qnum_select=qnum_select, set_ax=set_ax,
                                                   axspec=ax, figspec=fig, resolvingpower=resolvingpower, printout=printout,
                                                   color=color, visible=(not hide_comp),
                                                   lw=lw, removebase=removebase)
        if lam[0] >= lamRange[0] and lam[-1] <= lamRange[1]:
            d.append((lam,spec))
    if coadd and len(d) > 0:
        xy = coadd_spec(d)
        ax.plot(*xy, color=color_coadd, lw=lw_coadd)
        if returndata:
            return xy


def plot_spectra_dir_vel(files, vel_range=(-10,10), ySpecRange=(-0.1, 1e2), beamsize=None, printout=True,
    qnum_select=None, set_ax=True, figspecsize=(10,3)):
    vel_range = vel_range
    filename = files[0]
    h3 = fits.open(filename)
    ax, fig = plot_spectrum_combine(h3, xtype='vel', vel_range=vel_range, ySpecRange=ySpecRange, unit='Jy', printout=printout, set_ax=set_ax,
                                    xscale='linear', yscale='linear', visible=False,
                                    figspecsize=figspecsize, returnax=True, lw=1)

    d = []
    for filename in files:
        h3 = fits.open(filename)
        ax, fig = plot_spectrum_combine(h3, xtype='vel', returnax=True, printout=printout, qnum_select=qnum_select, set_ax=set_ax,
                                     axspec=ax, figspec=fig,
                                     lw=1, removebase=True, beamsize=beamsize)



def plot_contribution_function(c, nn=200, xmin=0, xmax=20, ymin=-10, ymax=10, nlev=100,
                               thr=1e-5, normalize=True,
                               return_data = False,
                               xscale='linear', yscale='linear',
                               figsize=(8,8)):
    from scipy.interpolate import griddata
    r = sqrt(c[:,0]**2 + c[:,1]**2)

    xi = linspace(xmin, xmax, nn)
    yi = linspace(ymin, ymax, nn)

    xi, yi = np.meshgrid(xi, yi)
    maxval = c[:,3].max()
    c[c[:,3] < maxval*thr, 3] = 0.0
    if normalize:
        c[:,3] = c[:,3] / maxval

    zi = griddata((r, c[:,2]), log10(c[:,3]), (xi, yi), method='linear')

    f = figure(figsize=figsize)
    pos = (0.15, 0.15, 0.8, 0.4)
    ax = f.add_axes(pos,
          xlabel='r (AU)',
          ylabel='z (AU)',
          autoscalex_on=False, autoscaley_on=False,
          xscale=xscale, yscale=yscale,
          xlim=(xmin, xmax), ylim=(ymin, ymax))

    C = ax.contourf(xi, yi, zi, nlev, cmap=mycm, extend='neither')
    for _ in C.collections:
        _.set_rasterized(True)
    #colorbar(C)
    set_axis_format(ax, graygrid=True, majorgridon=True, minorgridon=False)
    #
    idx = logical_not(isfinite(zi))
    zi[idx] = -99.0
    flx = sum(xi * 10**zi, axis=0)
    flx_t = sum(flx)
    flx = flx / flx_t
    fsum = cumsum(flx)
    pos = (0.15, 0.6, 0.8, 0.33)
    maxval = flx.max()
    maxtik = ceil(maxval/0.1)*0.1
    yticks = linspace(0,maxtik,num=int(maxtik/0.1)+1)
    ax = f.add_axes(pos,
          xlabel='', ylabel='Contribution',
          xticklabels=[],
          yticks=yticks,
          autoscalex_on=False, autoscaley_on=False,
          xscale=xscale, yscale=yscale,
          xlim=(xmin, xmax), ylim=(yticks[0], yticks[-1]))
    ax.plot(xi[0], flx,  lw=2, color='red')
    set_axis_format(ax, majorgridon=False, minorgridon=False)
    ax = twinx(ax)
    ax.plot(xi[0], fsum, lw=1, color='blue')
    ax.set_ylim((0,1))
    set_axis_format(ax, majorgridon=False, minorgridon=False)
    if return_data:
        return xi, yi, zi
    else:
        return

def plot_image(h, ext=2, data_in=None, edge_enhance=False,
               axmap=None, figmap=None, figmapsize=(6,6),
               xscale='linear', yscale='linear',
               xRange=None, yRange=None, logdata=False, valRange=None,
               xy_in_angle=False, beamsize=None):
    a = h[ext]
    n1 = a.header['NAXIS1']
    n2 = a.header['NAXIS2']
    dx = a.header['CDELT1']
    dy = a.header['CDELT2']
    k0x = a.header['CRPIX1']
    k0y = a.header['CRPIX2']
    x0 = a.header['CRVAL1']
    y0 = a.header['CRVAL2']
    
    x = (arange(1, n1+1) - k0x) * dx + x0
    y = (arange(1, n2+1) - k0y) * dy + y0
    
    if xy_in_angle:
        # In arcsec
        x /= h[0].header['DIST']
        y /= h[0].header['DIST']
        xyunit = '($^{\prime\prime}$)'
    else:
        xyunit = '(AU)'
    if data_in == None:
        data_in = h[ext].data
    if beamsize != None:
        from scipy.ndimage.filters import gaussian_filter as gsfl
        data = gsfl(data_in, ceil(beamsize*0.5/(dx/h[0].header['DIST'])))
    else:
        data = data_in
    if logdata:
        data = log10(data.copy())
    if edge_enhance:
        from scipy.ndimage.morphology import morphological_gradient
        data = morphological_gradient(data, size=(3,3))
    if figmap == None:
        figmap = figure(figsize=figmapsize)
    if axmap == None:
        if xRange == None:
            xRange = (x.min(), x.max())
            yRange = (y.min(), y.max())
        axmap = figmap.add_axes((0.1,0.1,0.85,0.85),
                              xlabel='X ' + xyunit,
                              ylabel='Y ' + xyunit,
                              autoscalex_on=False, autoscaley_on=False,
                              xscale=xscale, yscale=yscale,
                              xlim=xRange, ylim=yRange)
    #C = axmap. contourf(x, y, data, 100, cmap=cm.gray, )
    if valRange != None:
        vmin, vmax = valRange
    else:
        vmin, vmax = data.min(), data.max()
    C = axmap.imshow(data, interpolation='none', origin='lower', vmin=vmin, vmax=vmax,
                     extent=(x.min(), x.max(), y.min(), y.max()), cmap=cm.gray)
    colorbar(C)
    axmap.set_aspect('equal')
    set_axis_format(axmap, xscale=xscale, yscale=yscale,
                    majorgridon=False, minorgridon=False)
    
    
def extract_cube(h, dim0_range, dim0_free_range, data=None):
    if data == None:
        data = h[0].data
    return data[dim0_range[0]:dim0_range[1], :, :].mean(axis=0) - \
           data[dim0_free_range[0]:dim0_free_range[1], :, :].mean(axis=0)



def stepify(x, y):
    # Example: fill_between(*stepify(x, y))
    import numpy as np
    x, y = np.asarray(x), np.asarray(y)
    n = len(x)
    xnew = np.zeros(2*n)
    ynew = np.zeros(2*n)
    xnew[1:-1:2] = 0.5*(x[0:-1] + x[1:])
    xnew[2::2]   = 0.5*(x[0:-1] + x[1:])
    ynew[1:-1:2] = y[:-1]
    ynew[2::2]   = y[1:]
    xnew[0]  = x[0]  - 0.5*(x[1]-x[0])
    xnew[-1] = x[-1] + 0.5*(x[-1]-x[-2])
    ynew[0]  = y[0]
    ynew[-1] = y[-1]
    return xnew, ynew



def get_all_matches(el, word):
    m = []
    len_el = len(el)
    len_w = len(word)
    i = 0
    while(1):
        j = i + len_el
        if j > len_w:
            break
        if word[i:j] == el:
            m.append(i)
        i += 1
    return m

def get_leading_digits(word):
    s = ''
    for L in word:
        if '0' <= L and L <= '9':
            s += L
        else:
            break
    return s
    
def get_ele_counts(el, word, all_el):
    len_el = len(el)
    len_w = len(word)
    el_c = [c for c in all_el if (c.startswith(el) and c!=el)]
    idx = get_all_matches(el, word)
    ct = 0
    for i in idx:
        if len(filter(word[i:].startswith, el_c)) != 0:
            continue
        ii = i+len_el
        if i+len_el == len_w:
            ct += 1
            continue
        Ldig = get_leading_digits(word[ii:])
        if Ldig == '':
            ct += 1
            continue
        else:
            ct += int(Ldig)
    return ct

def get_molecule_names_from_file(fname, col_start=None):
    #col_start = 133
    with open(fname, 'r') as f:
        str_comment = f.readline()[1:].split()
    if col_start == None:
        col_start = str_comment.index('H')
    return str_comment[col_start:]

def update_elemental_abundance(ele, d, fname, all_elements):
    from numpy import zeros
    molecule_names = get_molecule_names_from_file(fname)
    n = len(d['H'])
    X_ele = zeros(n)
    for k in molecule_names:
        c = get_ele_counts(ele, k, all_elements)
        if c != 0:
            try:
                X_ele += c * d[k]
            except:
                raise KeyError(k)
    d['X['+ele+']'] = X_ele
    return

def update_density_of_X(X, d):
    d['n('+X+')'] = d[X] * d['n_gas']
    return


def update_stuff(ca, filename):
    update_density_of_X('H2O', ca)
    update_density_of_X('OH', ca)
    update_density_of_X('CO', ca)
    update_density_of_X('O', ca)
    update_elemental_abundance('O', ca, filename, all_elements)
    update_elemental_abundance('C', ca, filename, all_elements)
    ca['Tgas-Tdust'] = ca['Tgas'] - ca['Tdust']


def is_hydrocarbon(name):
    c_H = get_ele_counts('H', name, all_elements)
    c_C = get_ele_counts('C', name, all_elements)
    if c_H == 0 or c_C == 0:
        return False
    for el in all_elements:
        if el == 'H' or el == 'C':
            continue
        c_el = get_ele_counts(el, name, all_elements)
        if c_el > 0:
            return False
    return True


def has_nitrogen(name):
    c_N = get_ele_counts('N', name, all_elements)
    if c_N > 0:
        return True
    else:   
        return False



def get_hydrocarbons(all_names):
    hc = []
    for name in all_names:
        if is_hydrocarbon(name):
            hc.append(name)
    return hc


def get_nitrogens(all_names):
    ni = []
    for name in all_names:
        if has_nitrogen(name):
            ni.append(name)
    return ni



all_elements = ['H', 'D', 'He', 'C', 'N', 'O', 'Si', 'S', 'Fe', 'Na', 'Mg', 'Cl', 'P', 'F', 'Ne', 'Ar', 'K']


def query_by_approx_xy(key, x, y, d):
    idx = where(logical_and(logical_and(d['rmin'] <= x, x <= d['rmax']),
                            logical_and(d['zmin'] <= y, y <= d['zmax'])))[0]
    if len(idx) == 0:
        return None
    elif len(idx) == 1:
        return d[key][idx]
    else:
        raise KeyError('Cells overlap!')
