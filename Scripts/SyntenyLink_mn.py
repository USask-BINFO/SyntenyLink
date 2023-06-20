## Some sorting functions

### 1. To compare the whether two subgenomes are two different strands of the same genes

    #If true: Output is 1
    #If false: Output is 2

def compare_strands_subgenomes(string1, string2):
    # a function to test if two given subgenomes are forward or reverse strands
    if len(string1) == len(string2)+2 or len(string2)==len(string1)+2:
        count=0
        m=max(len(string1),len(string2))
        for i in range(m-2):
            # Compare characters at the same position
            if string1[i] == string2[i]:
                count=count+1
        if count==m-2:
            return 1
        else:
            return 0
    else:
        #they are not forward and reverse strand
        return 0
    # Iterate through three character in both strings

def compare_subgenomes(string1, string2):
  if string1==string2 or compare_strands_subgenomes(string1, string2):
     return 1
  else:
     return 0

### Function for Counting occureneces of a given chromosome in a given array
def count_occurrences(array, entry):
    if len(array) == 0:
        print("\n Empty array")
    count = 0
    for elem in array:
        if elem == entry:
            count += 1
    return count
# Example usage:
given_array = ["N1", "N2", "N3.r", "N4", "N2", "N5", "N2"]
given_entry = "N3.r"

occurrences = count_occurrences(given_array, given_entry)
# print("Number of occurrences:", occurrences)
    