#function to annotate a polygon id column in a dataframe with the polygon name

import sys
sys.path.insert(0, '/Users/david/github_datascience/vortexa/analytics/utils')
import sql_utils_david as sql_utils, math_utils




def annotate_polygons(df_in,column_df,df_poly='auto',column_poly=''):
 
 if (df_poly == 'auto'):
  pd.set_option("display.max_rows", 200)
  sql_engine = sql_utils.initiate_sql_engine(DB_NAME='backend_nov7')
  poly = pd.read_sql("""select * from new_polygons""")
 else:
  poly = df_poly
  
  
  
 if (column_poly == ''):
  cp = 'id'
 else:
  cp = column_poly
 annotated = df_in.join(poly.set_index(cp),on=column_df)
 #df_in.set_index(column_df).join(poly.set_index(column_poly))
 
 annotated_accepted = annotated[(annotated[cp].notnull())
 
 print(len(annotated) - len(annotated_accepted),' entries excluded from the output due to joining failure)
 print(100.*(len(annotated) - len(annotated_accepted))/len(annotated),' pc of total entries')
 
 return(annotated_accepted)
 
 

