import numpy as np
import pandas as pd




def annotate_polygons(df,polgons,searchstring):

 #identify polygon id corresponding to search string
 ref_id = polygons[polygons['id'].str.contains(searchstring)]['id']

 #return all polygons that contain the reference id in their parent id's
 #i.e all the children of the searchstring
 npoly = len(polygons)
 idx_in = []
 for idx in range(npoly):
  parent_id = polygons['parent_id'][idx]
  if (ref_id in parent_id):
   idx_in.append(idx)
   

 polygons_in = polygons[idx_in]
 
 
 #return all the df columns that contain any of the polygon id's in polygons_in
 id_in = df['id'].isin(polygons_in)

 return(df[id_in])