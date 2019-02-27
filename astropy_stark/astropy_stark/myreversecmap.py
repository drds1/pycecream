import matplotlib as mpl
def reverse_colourmap(cmap, name = 'my_cmap_r'):
     return mpl.colors.LinearSegmentedColormap(name, mpl.cm.revcmap(cmap._segmentdata))