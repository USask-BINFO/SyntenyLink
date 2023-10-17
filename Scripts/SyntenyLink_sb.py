import pandas as pd
import numpy as np

def update_df_synteny(C_df,synteny_file, df_subgenome_density_file_path, modified_synteny, main_bp_density_rmv):
    '''
    This function updates the df_synteny dataframe.

    Parameters
    ----------
    df_synteny : dataframe
        Dataframe of synteny.

    Returns
    -------
    df_synteny : dataframe
        Dataframe of synteny with chromosome names where the gene belongs.
    '''
        # Read the data
    df_subgenome_density = pd.read_excel(df_subgenome_density_file_path)

    # Get the column names from df_subgenome_density that start with N followed by a number
    column_names = df_subgenome_density.filter(regex=r'^N\d+').columns.tolist()

    # Replace the columns starting from the second column until the end of C_df_updated with the selected column names
    df_synteny = synteny_file.rename(columns={synteny_file.columns[i-1]: column_names[i-1] for i in range(1, len(column_names)+1)})

    # Start the index from 0
    # df_synteny.index = df_synteny.index + 1
    for i in range(len(main_bp_density_rmv)):
        for k in range(int(main_bp_density_rmv.iloc[i]['Row start #']), int(main_bp_density_rmv.iloc[i]['Row end #'] + 1)):
            for j in range(len(df_synteny.columns)):          
                if df_synteny.iloc[k, j] == 1 and df_synteny.columns[j] != main_bp_density_rmv.iloc[i]['subgenome1'] and df_synteny.columns[j] != main_bp_density_rmv.iloc[i]['subgenome2'] and df_synteny.columns[j] != main_bp_density_rmv.iloc[i]['subgenome3']:
                    df_synteny.iloc[k, j] = df_synteny.columns[j]
    #append the first column of the C_df to the df_synteny as first column
    df_synteny.insert(0, "Gene_id", C_df.iloc[:, 0])
    
    df_synteny.to_excel(modified_synteny)

    return df_synteny

def small_blk(df_focus, main_bp_density_rmv, chr_names_C):
    new_df_focus = pd.DataFrame(columns=['Row start #', 'Row end #', 'subgenome1', 'subgenome2', 'subgenome3'])
    new_index = 0
    sub1 = False
    sub2 = False
    sub3 = False

    for i in range(len(df_focus)):
        # only consider columns starting with N
        for j in range(2, len(df_focus.columns) - 1):
            # if df_focus.iloc[i]['Non_zero'] > 3:  
            if (df_focus.iloc[i][j] != 0 and df_focus.iloc[i][j] < 0.1):
                for s in range(len(main_bp_density_rmv)):
                    if main_bp_density_rmv.iloc[s]['Row start #'] == df_focus.iloc[i]['Row start #']:
                        if df_focus.iloc[i]['Non_zero'] - main_bp_density_rmv.iloc[s]['Non_zero'] < 2 :
                                # if column name j of df_focus is in the subgenome1 or subgenome2 or subgenome3 column of main_bp_density_rmv
                            for k in range(int(df_focus.iloc[i]['Row start #']), int(df_focus.iloc[i]['Row end #'] + 1)):
                                
                                new_df_focus.at[new_index, 'Row start #'] = int(df_focus.iloc[i]['Row start #'])
                                # take the row start where you first find C_new.iloc[k][j] == df_focus.columns[j]
                                if chr_names_C.iloc[k][j] == 1 or chr_names_C.iloc[k][j] == 0:
                                    continue

                                if chr_names_C.iloc[k][j] == df_focus.columns[j]:
                                    new_df_focus.at[new_index, 'Row end #'] = k-1
                                    new_df_focus.at[new_index, 'subgenome1'] = main_bp_density_rmv.iloc[s]['subgenome1']
                                    new_df_focus.at[new_index, 'subgenome2'] = main_bp_density_rmv.iloc[s]['subgenome2']
                                    new_df_focus.at[new_index, 'subgenome3'] = main_bp_density_rmv.iloc[s]['subgenome3']
                                    # add row start as the row start # of df_focus as a new row after i
                                    new_df_focus.loc[new_index + 1, 'Row start #'] = k
                                    # add row end as the row where it becomes not equal to df_focus.columns[j]
                                    row_end = k
                                    #window_size = 2
                                    while row_end < df_focus.loc[i]['Row end #']+1:
                                        if chr_names_C.iloc[row_end][j] == df_focus.columns[j]:
                                            chr_name = df_focus.columns[j]
                                            row_end += 1
                                        if chr_names_C.iloc[row_end][j] != df_focus.columns[j]:
                                            # chr_name = df_focus.columns[j]
                                            row_end += 1
                                            break
                                        break
                                        
                                    new_df_focus.at[new_index + 1, 'Row end #'] = row_end
                                    if chr_name.split('.')[0] == main_bp_density_rmv.iloc[s]['subgenome1'].split('.')[0]:
                                        if new_df_focus.at[new_index + 1, 'subgenome2'] != chr_name and new_df_focus.at[new_index + 1, 'subgenome3'] != chr_name:
                                            new_df_focus.at[new_index + 1, 'subgenome1'] = chr_name
                                            sub1 = True
                                    # if chr_name != main_bp_density_rmv.iloc[s]['subgenome1'] or chr_name != main_bp_density_rmv.iloc[s]['subgenome1'].split('.')[0]:
                                    #     new_df_focus.at[new_index + 1, 'subgenome1'] = main_bp_density_rmv.iloc[s]['subgenome1']
                                        # new_df_focus.at[new_index + 1, 'subgenome2'] = main_bp_density_rmv.iloc[s]['subgenome2']
                                        # new_df_focus.at[new_index + 1, 'subgenome3'] = main_bp_density_rmv.iloc[s]['subgenome3']
                                    if chr_name.split('.')[0] == main_bp_density_rmv.iloc[s]['subgenome2'].split('.')[0]:
                                        if new_df_focus.at[new_index + 1, 'subgenome1'] != chr_name and new_df_focus.at[new_index + 1, 'subgenome3'] != chr_name:
                                            new_df_focus.at[new_index + 1, 'subgenome2'] = chr_name
                                            sub2 = True

                                    if chr_name.split('.')[0] == main_bp_density_rmv.iloc[s]['subgenome3'].split('.')[0]:
                                        if new_df_focus.at[new_index + 1, 'subgenome1'] != chr_name and new_df_focus.at[new_index + 1, 'subgenome2'] != chr_name:
                                            new_df_focus.at[new_index + 1, 'subgenome3'] = chr_name
                                            sub3 = True
                                    if sub1 == False and sub2 == True and sub3 == True:
                                        new_df_focus.at[new_index + 1, 'subgenome1'] = main_bp_density_rmv.iloc[s]['subgenome1']
                                    if sub1 == True and sub2 == False and sub3 == True:
                                        new_df_focus.at[new_index + 1, 'subgenome2'] = main_bp_density_rmv.iloc[s]['subgenome2']
                                    if sub1 == True and sub2 == True and sub3 == False:
                                        new_df_focus.at[new_index + 1, 'subgenome3'] = main_bp_density_rmv.iloc[s]['subgenome3']
                                    if sub1 == False and sub2 == False and sub3 == True:
                                        new_df_focus.at[new_index + 1, 'subgenome1'] = main_bp_density_rmv.iloc[s]['subgenome1']
                                        new_df_focus.at[new_index + 1, 'subgenome2'] = main_bp_density_rmv.iloc[s]['subgenome2']
                                    if sub1 == False and sub2 == True and sub3 == False:
                                        new_df_focus.at[new_index + 1, 'subgenome1'] = main_bp_density_rmv.iloc[s]['subgenome1']
                                        new_df_focus.at[new_index + 1, 'subgenome3'] = main_bp_density_rmv.iloc[s]['subgenome3']
                                    if sub1 == True and sub2 == False and sub3 == False:
                                        new_df_focus.at[new_index + 1, 'subgenome2'] = main_bp_density_rmv.iloc[s]['subgenome2']
                                        new_df_focus.at[new_index + 1, 'subgenome3'] = main_bp_density_rmv.iloc[s]['subgenome3']
                                    if chr_name != main_bp_density_rmv.iloc[s]['subgenome1'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome2'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome3'] and chr_name != main_bp_density_rmv.iloc[s]['subgenome1'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome2'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome3'].split('.')[0]:
                                        new_df_focus.at[new_index + 1, 'subgenome1'] = chr_name
                                        if new_df_focus.at[new_index, 'subgenome2'] != new_df_focus.at[new_index + 1, 'subgenome1']:
                                            new_df_focus.at[new_index + 1, 'subgenome2'] = new_df_focus.at[new_index, 'subgenome2']
                                        if new_df_focus.at[new_index, 'subgenome3'] != new_df_focus.at[new_index + 1, 'subgenome1']:
                                            new_df_focus.at[new_index + 1, 'subgenome3'] = new_df_focus.at[new_index, 'subgenome3']
                                        if new_df_focus.at[new_index, 'subgenome2'] == new_df_focus.at[new_index + 1, 'subgenome1']:
                                            new_df_focus.at[new_index + 1, 'subgenome2'] = new_df_focus.at[new_index, 'subgenome1']
                                        if new_df_focus.at[new_index, 'subgenome3'] == new_df_focus.at[new_index + 1, 'subgenome1']:
                                            new_df_focus.at[new_index + 1, 'subgenome3'] = new_df_focus.at[new_index, 'subgenome1']
                                    if chr_name != main_bp_density_rmv.iloc[s]['subgenome2'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome1'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome3'] and chr_name != main_bp_density_rmv.iloc[s]['subgenome1'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome2'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome3'].split('.')[0]:
                                        new_df_focus.at[new_index + 1, 'subgenome2'] = chr_name
                                        if new_df_focus.at[new_index, 'subgenome1'] != new_df_focus.at[new_index + 1, 'subgenome2']:
                                            new_df_focus.at[new_index + 1, 'subgenome1'] = new_df_focus.at[new_index, 'subgenome1']
                                        if new_df_focus.at[new_index, 'subgenome3'] != new_df_focus.at[new_index + 1, 'subgenome2']:
                                            new_df_focus.at[new_index + 1, 'subgenome3'] = new_df_focus.at[new_index, 'subgenome3']
                                        if new_df_focus.at[new_index, 'subgenome1'] == new_df_focus.at[new_index + 1, 'subgenome2']:
                                            new_df_focus.at[new_index + 1, 'subgenome1'] = new_df_focus.at[new_index, 'subgenome2']
                                        if new_df_focus.at[new_index, 'subgenome3'] == new_df_focus.at[new_index + 1, 'subgenome2']:
                                            new_df_focus.at[new_index + 1, 'subgenome3'] = new_df_focus.at[new_index, 'subgenome2']
                                    if chr_name != main_bp_density_rmv.iloc[s]['subgenome3'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome1'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome2'] and chr_name != main_bp_density_rmv.iloc[s]['subgenome1'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome2'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome3'].split('.')[0]:
                                        new_df_focus.at[new_index + 1, 'subgenome3'] = chr_name
                                        if new_df_focus.at[new_index, 'subgenome1'] != new_df_focus.at[new_index + 1, 'subgenome3']:
                                            new_df_focus.at[new_index + 1, 'subgenome1'] = new_df_focus.at[new_index, 'subgenome1']
                                        if new_df_focus.at[new_index, 'subgenome2'] != new_df_focus.at[new_index + 1, 'subgenome3']:
                                            new_df_focus.at[new_index + 1, 'subgenome2'] = new_df_focus.at[new_index, 'subgenome2']
                                        if new_df_focus.at[new_index, 'subgenome1'] == new_df_focus.at[new_index + 1, 'subgenome3']:
                                            new_df_focus.at[new_index + 1, 'subgenome1'] = new_df_focus.at[new_index, 'subgenome3']
                                        if new_df_focus.at[new_index, 'subgenome2'] == new_df_focus.at[new_index + 1, 'subgenome3']:
                                            new_df_focus.at[new_index + 1, 'subgenome2'] = new_df_focus.at[new_index, 'subgenome3']
                                    if chr_name != main_bp_density_rmv.iloc[s]['subgenome1'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome1'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome2'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome3'] and chr_name != main_bp_density_rmv.iloc[s]['subgenome1'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome2'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome3'].split('.')[0]:
                                        new_df_focus.at[new_index + 1, 'subgenome1'] = chr_name
                                        if new_df_focus.at[new_index, 'subgenome2'] != new_df_focus.at[new_index + 1, 'subgenome1']:
                                            new_df_focus.at[new_index + 1, 'subgenome2'] = new_df_focus.at[new_index, 'subgenome2']
                                        if new_df_focus.at[new_index, 'subgenome3'] != new_df_focus.at[new_index + 1, 'subgenome1']:
                                            new_df_focus.at[new_index + 1, 'subgenome3'] = new_df_focus.at[new_index, 'subgenome3']
                                        if new_df_focus.at[new_index, 'subgenome2'] == new_df_focus.at[new_index + 1, 'subgenome1']:
                                            new_df_focus.at[new_index + 1, 'subgenome2'] = new_df_focus.at[new_index, 'subgenome1']
                                        if new_df_focus.at[new_index, 'subgenome3'] == new_df_focus.at[new_index + 1, 'subgenome1']:
                                            new_df_focus.at[new_index + 1, 'subgenome3'] = new_df_focus.at[new_index, 'subgenome1']
                                    if chr_name != main_bp_density_rmv.iloc[s]['subgenome2'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome1'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome2'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome3'] and chr_name != main_bp_density_rmv.iloc[s]['subgenome1'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome2'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome3'].split('.')[0]:
                                        new_df_focus.at[new_index + 1, 'subgenome2'] = chr_name
                                        if new_df_focus.at[new_index, 'subgenome1'] != new_df_focus.at[new_index + 1, 'subgenome2']:
                                            new_df_focus.at[new_index + 1, 'subgenome1'] = new_df_focus.at[new_index, 'subgenome1']
                                        if new_df_focus.at[new_index, 'subgenome3'] != new_df_focus.at[new_index + 1, 'subgenome2']:
                                            new_df_focus.at[new_index + 1, 'subgenome3'] = new_df_focus.at[new_index, 'subgenome3']
                                        if new_df_focus.at[new_index, 'subgenome1'] == new_df_focus.at[new_index + 1, 'subgenome2']:
                                            new_df_focus.at[new_index + 1, 'subgenome1'] = new_df_focus.at[new_index, 'subgenome2']
                                        if new_df_focus.at[new_index, 'subgenome3'] == new_df_focus.at[new_index + 1, 'subgenome2']:
                                            new_df_focus.at[new_index + 1, 'subgenome3'] = new_df_focus.at[new_index, 'subgenome2']
                                    if chr_name != main_bp_density_rmv.iloc[s]['subgenome3'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome1'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome2'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome3'] and chr_name != main_bp_density_rmv.iloc[s]['subgenome1'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome2'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome3'].split('.')[0]:
                                        new_df_focus.at[new_index + 1, 'subgenome3'] = chr_name
                                        if new_df_focus.at[new_index, 'subgenome1'] != new_df_focus.at[new_index + 1, 'subgenome3']:
                                            new_df_focus.at[new_index + 1, 'subgenome1'] = new_df_focus.at[new_index, 'subgenome1']
                                        if new_df_focus.at[new_index, 'subgenome2'] != new_df_focus.at[new_index + 1, 'subgenome3']:
                                            new_df_focus.at[new_index + 1, 'subgenome2'] = new_df_focus.at[new_index, 'subgenome2']
                                        if new_df_focus.at[new_index, 'subgenome1'] == new_df_focus.at[new_index + 1, 'subgenome3']:
                                            new_df_focus.at[new_index + 1, 'subgenome1'] = new_df_focus.at[new_index, 'subgenome3']
                                        if new_df_focus.at[new_index, 'subgenome2'] == new_df_focus.at[new_index + 1, 'subgenome3']:
                                            new_df_focus.at[new_index + 1, 'subgenome2'] = new_df_focus.at[new_index, 'subgenome3']
                                    if chr_name.split('.')[0] != main_bp_density_rmv.iloc[s]['subgenome1'].split('.')[0] and chr_name == None:
                                        new_df_focus.at[new_index + 1, 'subgenome1'] = new_df_focus.at[new_index, 'subgenome1']
                                    if chr_name.split('.')[0] != main_bp_density_rmv.iloc[s]['subgenome2'].split('.')[0] and chr_name == None:
                                        new_df_focus.at[new_index + 1, 'subgenome2'] = new_df_focus.at[new_index, 'subgenome2']
                                    if chr_name.split('.')[0] != main_bp_density_rmv.iloc[s]['subgenome3'].split('.')[0] and chr_name == None:
                                        new_df_focus.at[new_index + 1, 'subgenome3'] = new_df_focus.at[new_index, 'subgenome3']
                                    if row_end < df_focus.loc[i]['Row end #']:
                                        # add row start as the row where it becomes not equal to df_focus.columns[j]
                                        new_df_focus.at[new_index + 2, 'Row start #'] = row_end + 1
                                        new_df_focus.at[new_index + 2, 'Row end #'] = df_focus.loc[i]['Row end #']
                                        new_df_focus.at[new_index + 2, 'subgenome1'] = main_bp_density_rmv.iloc[s]['subgenome1']
                                        new_df_focus.at[new_index + 2, 'subgenome2'] = main_bp_density_rmv.iloc[s]['subgenome2']
                                        new_df_focus.at[new_index + 2, 'subgenome3'] = main_bp_density_rmv.iloc[s]['subgenome3']
                                        new_index += 3
                                        break
                                    break
    #                             if df_focus.iloc[i]['Row end #'] == len(chr_names_C):
    #                                 break
    #                                 break
                        if df_focus.iloc[i]['Non_zero'] - main_bp_density_rmv.iloc[s]['Non_zero'] >=2 :
                        # if column name j of df_focus is in the subgenome1 or subgenome2 or subgenome3 column of main_bp_density_rmv
                            for k in range(int(df_focus.iloc[i]['Row start #']), int(df_focus.iloc[i]['Row end #'] + 1)): 
                                new_df_focus.at[new_index, 'Row start #'] = int(df_focus.iloc[i]['Row start #'])                
                                # take the row start where you first find C_new.iloc[k][j] == df_focus.columns[j]
                                if chr_names_C.iloc[k][j] == 1 or chr_names_C.iloc[k][j] == 0:
                                    continue

                                if chr_names_C.iloc[k][j] == df_focus.columns[j]:
                                    
                                    new_df_focus.at[new_index, 'Row end #'] = k-1
                                    new_df_focus.at[new_index, 'subgenome1'] = main_bp_density_rmv.iloc[s]['subgenome1']
                                    new_df_focus.at[new_index, 'subgenome2'] = main_bp_density_rmv.iloc[s]['subgenome2']
                                    new_df_focus.at[new_index, 'subgenome3'] = main_bp_density_rmv.iloc[s]['subgenome3']
                                    # add row start as the row start # of df_focus as a new row after i
                                    new_df_focus.loc[new_index + 1, 'Row start #'] = k
                                    # add row end as the row where it becomes not equal to df_focus.columns[j]
                                    row_end = k
                                    #window_size = 2
                                    while row_end < df_focus.loc[i]['Row end #']+1:
                                        if chr_names_C.iloc[row_end][j] == df_focus.columns[j]:
                                            chr_name = df_focus.columns[j]
                                            row_end += 1
                                        if chr_names_C.iloc[row_end][j] != df_focus.columns[j]:
                                            # chr_name = df_focus.columns[j]
                                            row_end += 1
                                            break
                                        break
                                        
                                    new_df_focus.at[new_index + 1, 'Row end #'] = row_end
                                    if chr_name.split('.')[0] == main_bp_density_rmv.iloc[s]['subgenome1'].split('.')[0]:
                                        if new_df_focus.at[new_index + 1, 'subgenome2'] != chr_name and new_df_focus.at[new_index + 1, 'subgenome3'] != chr_name:
                                            new_df_focus.at[new_index + 1, 'subgenome1'] = chr_name
                                            sub1 = True
                                    # if chr_name != main_bp_density_rmv.iloc[s]['subgenome1'] or chr_name != main_bp_density_rmv.iloc[s]['subgenome1'].split('.')[0]:
                                    #     new_df_focus.at[new_index + 1, 'subgenome1'] = main_bp_density_rmv.iloc[s]['subgenome1']
                                        # new_df_focus.at[new_index + 1, 'subgenome2'] = main_bp_density_rmv.iloc[s]['subgenome2']
                                        # new_df_focus.at[new_index + 1, 'subgenome3'] = main_bp_density_rmv.iloc[s]['subgenome3']
                                    if chr_name.split('.')[0] == main_bp_density_rmv.iloc[s]['subgenome2'].split('.')[0]:
                                        if new_df_focus.at[new_index + 1, 'subgenome1'] != chr_name and new_df_focus.at[new_index + 1, 'subgenome3'] != chr_name:
                                            new_df_focus.at[new_index + 1, 'subgenome2'] = chr_name
                                            sub2 = True

                                    if chr_name.split('.')[0] == main_bp_density_rmv.iloc[s]['subgenome3'].split('.')[0]:
                                        if new_df_focus.at[new_index + 1, 'subgenome1'] != chr_name and new_df_focus.at[new_index + 1, 'subgenome2'] != chr_name:
                                            new_df_focus.at[new_index + 1, 'subgenome3'] = chr_name
                                            sub3 = True
                                    if sub1 == False and sub2 == True and sub3 == True:
                                        new_df_focus.at[new_index + 1, 'subgenome1'] = main_bp_density_rmv.iloc[s]['subgenome1']
                                    if sub1 == True and sub2 == False and sub3 == True:
                                        new_df_focus.at[new_index + 1, 'subgenome2'] = main_bp_density_rmv.iloc[s]['subgenome2']
                                    if sub1 == True and sub2 == True and sub3 == False:
                                        new_df_focus.at[new_index + 1, 'subgenome3'] = main_bp_density_rmv.iloc[s]['subgenome3']
                                    if sub1 == False and sub2 == False and sub3 == True:
                                        new_df_focus.at[new_index + 1, 'subgenome1'] = main_bp_density_rmv.iloc[s]['subgenome1']
                                        new_df_focus.at[new_index + 1, 'subgenome2'] = main_bp_density_rmv.iloc[s]['subgenome2']
                                    if sub1 == False and sub2 == True and sub3 == False:
                                        new_df_focus.at[new_index + 1, 'subgenome1'] = main_bp_density_rmv.iloc[s]['subgenome1']
                                        new_df_focus.at[new_index + 1, 'subgenome3'] = main_bp_density_rmv.iloc[s]['subgenome3']
                                    if sub1 == True and sub2 == False and sub3 == False:
                                        new_df_focus.at[new_index + 1, 'subgenome2'] = main_bp_density_rmv.iloc[s]['subgenome2']
                                        new_df_focus.at[new_index + 1, 'subgenome3'] = main_bp_density_rmv.iloc[s]['subgenome3']
                                    if chr_name != main_bp_density_rmv.iloc[s]['subgenome1'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome2'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome3'] and chr_name != main_bp_density_rmv.iloc[s]['subgenome1'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome2'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome3'].split('.')[0]:
                                        new_df_focus.at[new_index + 1, 'subgenome1'] = chr_name
                                        if new_df_focus.at[new_index, 'subgenome2'] != new_df_focus.at[new_index + 1, 'subgenome1']:
                                            new_df_focus.at[new_index + 1, 'subgenome2'] = new_df_focus.at[new_index, 'subgenome2']
                                        if new_df_focus.at[new_index, 'subgenome3'] != new_df_focus.at[new_index + 1, 'subgenome1']:
                                            new_df_focus.at[new_index + 1, 'subgenome3'] = new_df_focus.at[new_index, 'subgenome3']
                                        if new_df_focus.at[new_index, 'subgenome2'] == new_df_focus.at[new_index + 1, 'subgenome1']:
                                            new_df_focus.at[new_index + 1, 'subgenome2'] = new_df_focus.at[new_index, 'subgenome1']
                                        if new_df_focus.at[new_index, 'subgenome3'] == new_df_focus.at[new_index + 1, 'subgenome1']:
                                            new_df_focus.at[new_index + 1, 'subgenome3'] = new_df_focus.at[new_index, 'subgenome1']
                                    if chr_name != main_bp_density_rmv.iloc[s]['subgenome2'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome1'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome3'] and chr_name != main_bp_density_rmv.iloc[s]['subgenome1'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome2'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome3'].split('.')[0]:
                                        new_df_focus.at[new_index + 1, 'subgenome2'] = chr_name
                                        if new_df_focus.at[new_index, 'subgenome1'] != new_df_focus.at[new_index + 1, 'subgenome2']:
                                            new_df_focus.at[new_index + 1, 'subgenome1'] = new_df_focus.at[new_index, 'subgenome1']
                                        if new_df_focus.at[new_index, 'subgenome3'] != new_df_focus.at[new_index + 1, 'subgenome2']:
                                            new_df_focus.at[new_index + 1, 'subgenome3'] = new_df_focus.at[new_index, 'subgenome3']
                                        if new_df_focus.at[new_index, 'subgenome1'] == new_df_focus.at[new_index + 1, 'subgenome2']:
                                            new_df_focus.at[new_index + 1, 'subgenome1'] = new_df_focus.at[new_index, 'subgenome2']
                                        if new_df_focus.at[new_index, 'subgenome3'] == new_df_focus.at[new_index + 1, 'subgenome2']:
                                            new_df_focus.at[new_index + 1, 'subgenome3'] = new_df_focus.at[new_index, 'subgenome2']
                                    if chr_name != main_bp_density_rmv.iloc[s]['subgenome3'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome1'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome2'] and chr_name != main_bp_density_rmv.iloc[s]['subgenome1'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome2'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome3'].split('.')[0]:
                                        new_df_focus.at[new_index + 1, 'subgenome3'] = chr_name
                                        if new_df_focus.at[new_index, 'subgenome1'] != new_df_focus.at[new_index + 1, 'subgenome3']:
                                            new_df_focus.at[new_index + 1, 'subgenome1'] = new_df_focus.at[new_index, 'subgenome1']
                                        if new_df_focus.at[new_index, 'subgenome2'] != new_df_focus.at[new_index + 1, 'subgenome3']:
                                            new_df_focus.at[new_index + 1, 'subgenome2'] = new_df_focus.at[new_index, 'subgenome2']
                                        if new_df_focus.at[new_index, 'subgenome1'] == new_df_focus.at[new_index + 1, 'subgenome3']:
                                            new_df_focus.at[new_index + 1, 'subgenome1'] = new_df_focus.at[new_index, 'subgenome3']
                                        if new_df_focus.at[new_index, 'subgenome2'] == new_df_focus.at[new_index + 1, 'subgenome3']:
                                            new_df_focus.at[new_index + 1, 'subgenome2'] = new_df_focus.at[new_index, 'subgenome3']
                                    if chr_name != main_bp_density_rmv.iloc[s]['subgenome1'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome1'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome2'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome3'] and chr_name != main_bp_density_rmv.iloc[s]['subgenome1'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome2'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome3'].split('.')[0]:
                                        new_df_focus.at[new_index + 1, 'subgenome1'] = chr_name
                                        if new_df_focus.at[new_index, 'subgenome2'] != new_df_focus.at[new_index + 1, 'subgenome1']:
                                            new_df_focus.at[new_index + 1, 'subgenome2'] = new_df_focus.at[new_index, 'subgenome2']
                                        if new_df_focus.at[new_index, 'subgenome3'] != new_df_focus.at[new_index + 1, 'subgenome1']:
                                            new_df_focus.at[new_index + 1, 'subgenome3'] = new_df_focus.at[new_index, 'subgenome3']
                                        if new_df_focus.at[new_index, 'subgenome2'] == new_df_focus.at[new_index + 1, 'subgenome1']:
                                            new_df_focus.at[new_index + 1, 'subgenome2'] = new_df_focus.at[new_index, 'subgenome1']
                                        if new_df_focus.at[new_index, 'subgenome3'] == new_df_focus.at[new_index + 1, 'subgenome1']:
                                            new_df_focus.at[new_index + 1, 'subgenome3'] = new_df_focus.at[new_index, 'subgenome1']
                                    if chr_name != main_bp_density_rmv.iloc[s]['subgenome2'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome1'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome2'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome3'] and chr_name != main_bp_density_rmv.iloc[s]['subgenome1'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome2'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome3'].split('.')[0]:
                                        new_df_focus.at[new_index + 1, 'subgenome2'] = chr_name
                                        if new_df_focus.at[new_index, 'subgenome1'] != new_df_focus.at[new_index + 1, 'subgenome2']:
                                            new_df_focus.at[new_index + 1, 'subgenome1'] = new_df_focus.at[new_index, 'subgenome1']
                                        if new_df_focus.at[new_index, 'subgenome3'] != new_df_focus.at[new_index + 1, 'subgenome2']:
                                            new_df_focus.at[new_index + 1, 'subgenome3'] = new_df_focus.at[new_index, 'subgenome3']
                                        if new_df_focus.at[new_index, 'subgenome1'] == new_df_focus.at[new_index + 1, 'subgenome2']:
                                            new_df_focus.at[new_index + 1, 'subgenome1'] = new_df_focus.at[new_index, 'subgenome2']
                                        if new_df_focus.at[new_index, 'subgenome3'] == new_df_focus.at[new_index + 1, 'subgenome2']:
                                            new_df_focus.at[new_index + 1, 'subgenome3'] = new_df_focus.at[new_index, 'subgenome2']
                                    if chr_name != main_bp_density_rmv.iloc[s]['subgenome3'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome1'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome2'] and chr_name != new_df_focus.at[new_index + 1, 'subgenome3'] and chr_name != main_bp_density_rmv.iloc[s]['subgenome1'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome2'].split('.')[0] and chr_name != main_bp_density_rmv.iloc[s]['subgenome3'].split('.')[0]:
                                        new_df_focus.at[new_index + 1, 'subgenome3'] = chr_name
                                        if new_df_focus.at[new_index, 'subgenome1'] != new_df_focus.at[new_index + 1, 'subgenome3']:
                                            new_df_focus.at[new_index + 1, 'subgenome1'] = new_df_focus.at[new_index, 'subgenome1']
                                        if new_df_focus.at[new_index, 'subgenome2'] != new_df_focus.at[new_index + 1, 'subgenome3']:
                                            new_df_focus.at[new_index + 1, 'subgenome2'] = new_df_focus.at[new_index, 'subgenome2']
                                        if new_df_focus.at[new_index, 'subgenome1'] == new_df_focus.at[new_index + 1, 'subgenome3']:
                                            new_df_focus.at[new_index + 1, 'subgenome1'] = new_df_focus.at[new_index, 'subgenome3']
                                        if new_df_focus.at[new_index, 'subgenome2'] == new_df_focus.at[new_index + 1, 'subgenome3']:
                                            new_df_focus.at[new_index + 1, 'subgenome2'] = new_df_focus.at[new_index, 'subgenome3']
                                    if chr_name.split('.')[0] != main_bp_density_rmv.iloc[s]['subgenome1'].split('.')[0] and chr_name == None:
                                        new_df_focus.at[new_index + 1, 'subgenome1'] = new_df_focus.at[new_index, 'subgenome1']
                                    if chr_name.split('.')[0] != main_bp_density_rmv.iloc[s]['subgenome2'].split('.')[0] and chr_name == None:
                                        new_df_focus.at[new_index + 1, 'subgenome2'] = new_df_focus.at[new_index, 'subgenome2']
                                    if chr_name.split('.')[0] != main_bp_density_rmv.iloc[s]['subgenome3'].split('.')[0] and chr_name == None:
                                        new_df_focus.at[new_index + 1, 'subgenome3'] = new_df_focus.at[new_index, 'subgenome3']
                                    if row_end < df_focus.loc[i]['Row end #']:
                                        # add row start as the row where it becomes not equal to df_focus.columns[j]
                                        new_df_focus.at[new_index + 2, 'Row start #'] = row_end + 1
                                        new_df_focus.at[new_index + 2, 'Row end #'] = df_focus.loc[i]['Row end #']
                                        new_df_focus.at[new_index + 2, 'subgenome1'] = main_bp_density_rmv.iloc[s]['subgenome1']
                                        new_df_focus.at[new_index + 2, 'subgenome2'] = main_bp_density_rmv.iloc[s]['subgenome2']
                                        new_df_focus.at[new_index + 2, 'subgenome3'] = main_bp_density_rmv.iloc[s]['subgenome3']
                                        new_index += 3
                                        break
                                    break
    return new_df_focus
                
                    

def multi_col(new_df_focus):
    for l in range(len(new_df_focus)):
        for k in range(l+1, len(new_df_focus)):
            
            if new_df_focus.iloc[l]['Row start #'] == new_df_focus.iloc[k]['Row start #']:
                new_df_focus.iloc[k]['Row start #'] = new_df_focus.iloc[k-1]['Row end #']+1

    return new_df_focus

def all_blk(new_df_focus, main_bp_density_rmv):

    new_sbp = pd.DataFrame(columns=['Row start #', 'Row end #', 'subgenome1', 'subgenome2', 'subgenome3'])
    l = 0

    # Add entries from main_bp_density_rmv
    for i in range(len(new_df_focus)):
        if i == 0 and new_df_focus.iloc[i]['Row start #'] != 0:
            entry = new_df_focus.iloc[i]['Row start #']
            # get the index of the row where the row start # is equal to the entry in main_bp_density_rmv and appen all the rows before that to new_sbp
            index = main_bp_density_rmv[main_bp_density_rmv['Row start #'] == entry].index.values.astype(int)[0]
            for j in range(index):
                new_sbp.at[l, 'Row start #'] = main_bp_density_rmv.iloc[j]['Row start #']
                new_sbp.at[l, 'Row end #'] = main_bp_density_rmv.iloc[j]['Row end #']
                new_sbp.at[l, 'subgenome1'] = main_bp_density_rmv.iloc[j]['subgenome1']
                new_sbp.at[l, 'subgenome2'] = main_bp_density_rmv.iloc[j]['subgenome2']
                new_sbp.at[l, 'subgenome3'] = main_bp_density_rmv.iloc[j]['subgenome3']
                l += 1

        if i == 0 and new_sbp.iloc[0]['Row start #'] == 0:
            k = l
        if i < len(new_df_focus)-1:
            if new_df_focus.iloc[i]['Row end #'] + 1 == new_df_focus.iloc[i + 1]['Row start #']:
                    new_sbp.at[k, 'Row start #'] = new_df_focus.iloc[i]['Row start #']
                    new_sbp.at[k, 'Row end #'] = new_df_focus.iloc[i]['Row end #']
                    new_sbp.at[k, 'subgenome1'] = new_df_focus.iloc[i]['subgenome1']
                    new_sbp.at[k, 'subgenome2'] = new_df_focus.iloc[i]['subgenome2']
                    new_sbp.at[k, 'subgenome3'] = new_df_focus.iloc[i]['subgenome3']
                    k += 1

            if new_df_focus.iloc[i]['Row end #'] + 1 != new_df_focus.iloc[i + 1]['Row start #']:
                new_sbp.at[k, 'Row start #'] = new_df_focus.iloc[i]['Row start #']
                new_sbp.at[k, 'Row end #'] = new_df_focus.iloc[i]['Row end #']
                new_sbp.at[k, 'subgenome1'] = new_df_focus.iloc[i]['subgenome1']
                new_sbp.at[k, 'subgenome2'] = new_df_focus.iloc[i]['subgenome2']
                new_sbp.at[k, 'subgenome3'] = new_df_focus.iloc[i]['subgenome3']
                k += 1
                for s in range(len(main_bp_density_rmv)):

                    if new_sbp.iloc[k - 1]['Row end #'] + 1 == main_bp_density_rmv.iloc[s]['Row start #'] and main_bp_density_rmv.iloc[s]['Row start #'] != new_df_focus.iloc[i + 1]['Row start #']:
                        new_sbp.at[k, 'Row start #'] = main_bp_density_rmv.iloc[s]['Row start #']
                        new_sbp.at[k, 'Row end #'] = main_bp_density_rmv.iloc[s]['Row end #']
                        new_sbp.at[k, 'subgenome1'] = main_bp_density_rmv.iloc[s]['subgenome1']
                        new_sbp.at[k, 'subgenome2'] = main_bp_density_rmv.iloc[s]['subgenome2']
                        new_sbp.at[k, 'subgenome3'] = main_bp_density_rmv.iloc[s]['subgenome3']
                        k += 1
                        # break
                    if main_bp_density_rmv.iloc[s]['Row end #'] + 1 == new_df_focus.iloc[i + 1]['Row start #']:
                        # new_sbp.at[k, 'Row start #'] = new_df_focus.iloc[i+1]['Row start #'] 
                        # new_sbp.at[k, 'Row end #'] = new_df_focus.iloc[i+1]['Row end #'] 
                        # new_sbp.at[k, 'subgenome1'] = new_df_focus.iloc[i+1]['subgenome1']
                        # new_sbp.at[k, 'subgenome2'] = new_df_focus.iloc[i+1]['subgenome2']
                        # new_sbp.at[k, 'subgenome3'] = new_df_focus.iloc[i+1]['subgenome3']
                        # k += 1
                        break
                # break
        if i == len(new_df_focus)-1:
            new_sbp.at[k, 'Row start #'] = new_df_focus.iloc[i]['Row start #']
            new_sbp.at[k, 'Row end #'] = new_df_focus.iloc[i]['Row end #']
            new_sbp.at[k, 'subgenome1'] = new_df_focus.iloc[i]['subgenome1']
            new_sbp.at[k, 'subgenome2'] = new_df_focus.iloc[i]['subgenome2']
            new_sbp.at[k, 'subgenome3'] = new_df_focus.iloc[i]['subgenome3']
            k += 1
            if new_sbp.at[k-1, 'Row end #'] +1 == main_bp_density_rmv.iloc[-1]['Row start #']:
                new_sbp.at[k, 'Row start #'] = main_bp_density_rmv.iloc[-1]['Row start #']
                new_sbp.at[k, 'Row end #'] = main_bp_density_rmv.iloc[-1]['Row end #']
                new_sbp.at[k, 'subgenome1'] = main_bp_density_rmv.iloc[-1]['subgenome1']
                new_sbp.at[k, 'subgenome2'] = main_bp_density_rmv.iloc[-1]['subgenome2']
                new_sbp.at[k, 'subgenome3'] = main_bp_density_rmv.iloc[-1]['subgenome3']
                k += 1
            elif new_sbp.at[k-1, 'Row end #'] == main_bp_density_rmv.iloc[-1]['Row start #']:
                new_sbp.at[k, 'Row start #'] = main_bp_density_rmv.iloc[-1]['Row start #']+1
                new_sbp.at[k, 'Row end #'] = main_bp_density_rmv.iloc[-1]['Row end #']
                new_sbp.at[k, 'subgenome1'] = main_bp_density_rmv.iloc[-1]['subgenome1']
                new_sbp.at[k, 'subgenome2'] = main_bp_density_rmv.iloc[-1]['subgenome2']
                new_sbp.at[k, 'subgenome3'] = main_bp_density_rmv.iloc[-1]['subgenome3']
                k += 1
            elif new_sbp.at[k-1, 'Row end #'] != main_bp_density_rmv.iloc[-1]['Row start #']:
                new_sbp.at[k, 'Row start #'] = new_sbp.at[k-1, 'Row end #']+1
                new_sbp.at[k, 'Row end #'] = main_bp_density_rmv.iloc[-1]['Row end #']
                new_sbp.at[k, 'subgenome1'] = main_bp_density_rmv.iloc[-1]['subgenome1']
                new_sbp.at[k, 'subgenome2'] = main_bp_density_rmv.iloc[-1]['subgenome2']
                new_sbp.at[k, 'subgenome3'] = main_bp_density_rmv.iloc[-1]['subgenome3']
                k += 1
    return new_sbp
            




