# FromGRNstoCTTMs

The central top-level functions for these analyses appear in the libraries Organismo2 and Spaces. These two libraries as well as others with auxiliary functions are in the libs directory.  get\_att.cc and random\_walk.cc are two example programs that call functions in these libraries.

## Organismo2
### Main Functions:
```text
Organismo2
├── make_org() -- It gets all attractors of an N-gene GRN and constructs its CTTM. 
│   ├── initial_vector() -- It sets the initial activity state, in which all genes are turned off.
│   ├── contador() -- It reserves memory space.
│   ├── num_hileras() -- It reserves memory space.
│   ├── new_find_all_att() -- It gets all the cell types from a GRN.
│   └── new_fill_diff_cell_nw() -- It makes the CTTM with the cell types and the possible transitions. 
├── make_multi_org() -- It gets all attractors of multiple N-gene GRNs and constructs their CTTMs. It saves  
│   │                   all the cell types in a txt file. It deletes GRN and CTTM every iteration.
│   ├── initial_vector() --  It sets the initial activity state, in which all genes are turned off.
│   ├── contador() -- It reserves memory space.
│   ├── num_hileras() -- It reserves memory space.
│   ├── new_find_all_att() -- It gets all the cell types from a GRN.
│   ├── new_fill_diff_cell_nw() -- It makes the CTTM with the cell types and the possible transitions. 
│   ├── delete_grn() -- It deletes GRN.
│   ├── delete_diff() -- It deletes CTTM.
│   └── get_all_att() -- It saves all cell types in a txt file.
├── make_org_zero() -- It gets all attractors of an N-gene GRN and constructs its CTTM. 
│   │                  It checks whether the downstream cell type has all its genes inactive.
│   ├── initial_vector() --  It sets the initial activity state, in which all genes are turned off.
│   ├── contador() -- It reserves memory space.
│   ├── num_hileras() -- It reserves memory space.
│   ├── new_find_all_att() -- It gets all the cell types from a GRN.
│   ├── new_fill_diff_cell_nw() -- It makes the CTTM with the cell types and the possible transitions. 
│   ├── zero_att_position() -- It defines a CTTM target in which the downstream cell type has all its genes inactive and
│   │                          compares the values between the CTTM and the target.
│   ├── delete_grn() -- It deletes GRN.
│   └──  delete_diff() -- It deletes CTTM.
```


## Spaces
### Main Functions:
```text
├── make_all_nw() -- It makes all the GRNs with N genes, excludes isomorphisms and
│   │                saves their properties.
│   ├── get_n_nw()
│   ├── func_combination()
│   ├── return_grn()
│   ├── return_matriz_nodo()
│   ├── return_nw_positive_interactions()
│   ├── check_interaction_rows()
│   ├── check_nodes_interaction()
│   ├── check_one_component()
│   ├── get_positive_interactions()
│   ├── get_nwgrn_db()
│   └── nw_combination()
├── get_distance_db() -- It gets the number of changes required to switch from one network to another.
│   │ 
│   ├── func_combination()
│   ├── get_nwgrn_original()
│   └── nw_combination_dif()
├── comp_all_space() -- It tests all posible GRNs with N genes and compares the CTTM with the function's
│   │                   CTTM target.
│   ├── get_n_nw()
│   ├── func_combination()
│   ├── nw_combination_dif()
│   └── save_grn()
├── comp_two_db() -- It compares two GRN data bases and saves in a third data base all non-isomorphic GRN.
│   │                
│   ├── func_combination() 
│   ├── del_1row()
│   ├── get_nwgrn_original()
│   ├── return_grn()
│   ├── get_nwgrn_original()
│   ├── get_nwgrn_db()
│   └── nw_combination()
├── del_1row() -- It deletes the first network in a data base if it has mora than one SCC or if
│   │             the number of interactions per row does not decrease.
│   ├── get_nwgrn_original()
│   ├── check_one_component()
│   └── check_interaction_rows()
├── red_db() -- It reduces the size of a data base dividing it into equal parts.
│   │ 
│   ├── func_combination()
│   ├── get_nwgrn_original()
│   ├── return_grn()
│   ├── get_nwgrn_original()
│   ├── get_nwgrn_db()
│   ├── nw_combination()
│   └── return_grn()
├── make_random_wells() -- The equivalent of a random walk in GRNs in which every step is
│   │                       a change in the interaction.
│   └── save_grn()
├── space_wells_txt() -- It makes groups of GRNs checking which of them are just a step away in a data base (txt)
│   │                    using a flag to identify them.
│   ├── fill_false()
│   ├── get_nwgrn_original()
│   ├── change_state()
│   ├── random_mutation()
│   └── check_mutation() 
├── space_wells() -- It makes groups of GRNs checking which of them are just a step away in a data base.
│   │                    
│   ├── get_nwgrn_original()
│   ├── get_nwgrn_memory()
│   ├── random_mutation()
│   ├── review_mutation()
│   └──  change_all_well()
├── random_walks() -- Beginning with a GRN that makes a CTTM target, it makes a random change. If the CTTM 
│   │                 is the target we move to that GRN.
│   ├── get_nwgrn_original()
│   ├── random_mutation()
│   └── save_grn()
├── random_walks_2ears() -- Similar to random_walks() but there is more than one CTTM target. Therefore, it can walk 
│   │                       traverse different CTTMs.
│   ├── get_nwgrn_original()
│   ├── random_mutation()
│   └── save_grn()
├── random_nw() -- It makes random GRNs and compares the CTTM with a CTTM target.
│   │ 
│   └── save_grn()
├── review_mutation() -- It gets grn_check and compares with grn to see if they're the same.
│   │
│   └── get_nwgrncheck_memory()
├──  check_mutation() -- It checks if a GRN are in a data base.
│   │ 
│   └── get_nwgrn_db()
├── count_one_step() -- Prints the network index a step away.
│   │ 
│   ├── get_nwgrn_original()
│   ├── get_nwgrn_db()
│   └── number_of_diff() 
├──  count_number_differences() -- It compares a GRN with all GRNs in a data base and checks the minimum distance  
│   │                              (number of changes) between them.
│   ├── get_nwgrn_original()
│   └── get_nwgrn_db() 
├──  get_1fb_clusters() -- It returns the number of positive and negative auto interactions in a GRN
│   │ 
│   └── get_nwgrn_original()
├──  get_2fb_clusters() -- It seeks positive feedback loops and gets the auto interactions and interactions with
│   │                      other genes.
│   ├── get_nwgrn_original() 
│   ├── get_numb_conec_extra()
│   └── get_auto_int()
├──  get_3fb_clusters() -- It returns the number of 3 genes positive feedback loops in a GRN.
│   │ 
│   └── get_nwgrn_original()
├──  get_4fb_clusters() -- It returns the number of 4 genes positive feedback loops in a GRN.
│   │ 
│   └── get_nwgrn_original()
├──  remove_cicatt() -- It removes all CTTM whith SCC.
│   │                     
│   └── get_nwgrn_original()
├──  genes_comp() -- It gets the number of interactions between genes in the same feedback loop.
│   │ 
│   ├── get_group()
│   └── get_nwgrn_original()

```
### Helper functions:
```text

    return_grn() -- It returns one value from the GRN matrix. 
    return_matriz_nodo() -- It returns the number of interactions.
    return_nw_positive_interactions() -- It returns the number of positive interactions.
    get_n_nw() -- If GRN is in a file as a column, it saves it in GraphI0.
    get_positive_interactions() -- It takes a line from a file and assigns it to num_positive_interactions_db.
    get_nwgrn_db() -- It assigns a GRN to grn_check.
    get_nwgrn_original() -- It assigns a GRN to grn.
    save_grn() -- It saves a GRN in a txt file.
    check_one_component() -- It checks if a GRN has more than one SCC.
    check_interaction_rows() -- It checks if the number of interactions in a matrix decreases every row.
    check_nodes_interaction() -- It checks if there are no empty rows.
    nw_combination() -- It makes the GRN combinations to check isomorphisms. 
    func_combination() -- It makes the list of combinations showing how to move the genes.
    check_nwdb_interactions() -- It counts the number of interactions from a saved network.
    review_tar_nw() -- It counts the number of positive feedback loops between genes in GRNs.
    nw_combination_dif() -- It makes the CTTM combinations to check ismorphisms.
    change_all_well() -- It changes the flags in a GRN+flag matrix for other flag.
    get_nwgrn_memory() -- It makes a GRN from a matrix row.
    get_nwgrncheck_memory() -- It makes a GRN_check from a matrix row.
    fill_false() -- It copies a data base and add a flag at the end of every line.
    change_state() -- It changes the flag 0 to 1 at the end of every line.
    random_mutation() -- It changes randomly a matrix value and save the original value.
    get_next_idx() -- It finds the next GRN which not belong to a gruop using a flag.
    number_of_diff() -- It counts the number of different entries between grn and grn_check.
    get_numb_conec_extra() -- It gets the interactions between genes in a feedback loop and other genes.
    get_auto_int() -- It gets the interactions between genes in a same feedback loop.
    get_list_1scc() -- It gets the index from all GRN with N genes and one SCC.
    get_2fb_members() -- List of all genes that are part of a 2 genes feedback loop.
    get_3fb_members() -- List of all genes that are part of a 3 genes feedback loop.
    get_4fb_members() -- List of all genes that are part of a 4 genes feedback loop.
    conection_gen_extra() -- Interactions between genes that belong to a feedback loop as opposed to those that do.
    node_negative_connections_only() -- It seeks if a gene does not belong to a feedback loop get negative interactions.
    all_nodes_are_fb() -- It seeks if all genes in a GRN belong to a feedback loop.
    two_three_nodes() -- It seeks for feedback loops of 2,3 y 4 genes size and interactions with genes dont belong.
    rev_nodes() -- It seeks the number of genes with positive interactions and the sum of values.
    space_distribution() -- It counts the number of differents in every network in a data base.
    seek_cicatt_allspace() -- It seeks for cyclic attractors.
    make_graphviz() -- Print all the neworks in a data base.
    get_group() -- It gets the genes in the same feedback loop from the cell types. 
```

