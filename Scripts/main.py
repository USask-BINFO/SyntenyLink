# Description: This is the main file for the SyntenyLink tool. It will call the main function from the accuracy_ADT.py file and run the main function.

import accuracy_ADT as A



# Create a breakpoint class object
synteny = A.Acc()

# Run the main function
subgenomes_final, GT, n_sub, swapped_sub, subgenome_letters = synteny.main()
# Get the final subgenome placements
comb_all = synteny.all_possible_combinations(subgenomes_final, n_sub, swapped_sub, GT, subgenome_letters)

# if ground truth is available
if GT is not None:
    # Get the still missing genes
    synteny.find_missing_genes(GT, subgenomes_final)
