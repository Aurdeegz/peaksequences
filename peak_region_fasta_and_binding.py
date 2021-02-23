"""
Finding the binding sites
"""


########################################################################
#
#          Imports

import string   # Module that helps with string methods. Used here
                # for replacing a substring with the capitalized
                # version of that substring.
import os       # Used for checking if a filepath exists

#
#
########################################################################
#
#    Pre-defined variables

# This dictionary contains keys that are file extensions
# and values of dictionaries with information on which
# column of the data file will have the desired information.
                   # narrowPeak files (from MACS3):
peakfile_formats ={"narrowPeak" : {"chrom" : "1",  # Chromosome, column1
                                   "start" : "2",  # Region start, column2
                                   "end"   : "3"}, # Region end, column3
                   # Excel spreadsheet files (from MACS3)
                   "xls"        : {"chrom" : "1",  # Chromosome, column1
                                   "start" : "2",  # Region start, column2
                                   "end"   : "3"}} # Region end, column3


# This is a list of strings containing the single letter identifiers
# for each nucleotide. "N" is the symbol for an unidentified nucleotide
# in a sequence.
nukes = ["A", "T", "C", "G", "N"]

#
#
########################################################################
#
#    Helper Functions

def str_a_list(a_list):
    """
    Given a list, return all elements of the list
    as a single string.
    """
    # Assure that the user gave a list object
    assert type(a_list) == list, "The input was not a list."
    # Initialize the output string
    outstr = ""
    # Loop over the elements in the list
    for element in a_list:
        # Update the outstr with the new element at the end
        outstr = f"{outstr}{element}"
    # Return the outstr variable
    return outstr

def reverse_a_str(a_string):
    """
    Given a string, return that string with the
    characters in reverse order
    """
    # Assure that the user input a string
    assert type(a_string) == str, "The input was not a string"
    # Turn the string into a list. Each character will be
    # a unique element of the list
    str_list = list(a_string)
    # Use list comprehension to iterate over the elements
    # of the list in reverse order. 
    # Note: a_list[-i] gets the ith element from the end of the list. 
    str_list = [str_list[-i] for i in range(1,len(str_list) + 1)]
    # Use the str_a_list() function to turn the list into a
    # string and return that string.
    return str_a_list(str_list)

def add_newline_every_x_chars(string, x = 80):
    """
    Given a string and an integer number x (default is 80),
    return the string with a newline character every x 
    """
    # Assuer that the user input a string and an integer
    assert type(string) == str, f"The input {string} is not a string"
    assert type(x) == int, f"The input {x} is not an integer type"
    # Assure that the value of x is less than or equal to
    # the length of the string
    assert len(string) >= x, f"The distance {x} is greater than the size of the string"
    # Initialize the newstring variable as an empty string
    newstring = ""
    # Initialize the count variable as 0
    count=0
    # Loop over the characters in the string
    for char in string:
        # If the count is zero
        if count == 0:
            # Then update the newstring with this character
            newstring=f"{newstring}{char}"
            # and increase the count
            count+=1
        # Or if the count is not a multiple of the x integer
        # (count % x is the modular operator and returns the
        # remainder of dividing count by x)
        elif count % x != 0:
            # Then update the newstring variable
            newstring=f"{newstring}{char}"
            # and increase the count
            count+=1
        # Otherwise, the count is a multiple of x. Thus,
        else:
            # add the current character and a newline character
            # to the newstring
            newstring=f"{newstring}{char}\n"
            # and reset the count to zero
            count=0
    # At the end, return the newstring variable.
    return newstring

#
#
#########################################################################
#
#      Frequency matrix calculator: Usage is optional

def check_motif_lengths(motifs_list):
    """
    Given a list of (motif, strand) tuples, check the length of the motifs
    and partition them based on their lengths
    """
    # Initialize the dictionary to hold partitions
    motifs_dict = {}
    # Loop over the motifs in the motif list
    for motif in motifs_list:
        # Check to see if the motif length is in the motifs_dict key
        # values already. If not
        if f"{len(motif[0])}" not in motifs_dict.keys():
            # Then initialize the motifs dictionary with this key
            # and a list as the value
            motifs_dict[f'{len(motif[0])}'] = [motif]
        # If the length of the motif is already a key in the dictionary
        else:
            # Then simply update the list with this sequence
            motifs_dict[f'{len(motif[0])}'].append(motif)
    # Once the loop is completed, return the motifs dictionary.
    return motifs_dict

def get_frequency_dict(motifs_list):
    """
    Given a list of (motif, strand) tuples, return a dictioanry
    counting the frequency of each nucleotide in each position
    of the sequence
    """
    assert type(motifs_list) == list, "The motifs should be a list of tuples..."
    for motif in motifs_list:
        assert type(motif) == tuple, "The motifs should be a list of tuples..."
        assert type(motif[0]) == str, "The zeroeth element of each tyuple should be a str type"
        assert motif[1] == "+" or motif[1] == "-", "The first element of each tuple should denote strandedness: + or -"
    # Initialize the frequency_dict variable with the keys as
    # nucleotides and empty dictionaries as values.
    frequency_dict = {"A": {},
                      "G": {},
                      "T": {},
                      "C": {}}
    # Loop over the motifs in the motifs_list
    for motif in motifs_list:
        # Loop over the siz of the motifs
        for i in range(len(motif[0])):
            # Loop over the keys and values of the frequency dictionary
            for key, value in frequency_dict.items():
                # If the position is not already a key in the subdictionary,
                if f"{i}" not in value.keys():
                    # and if the current nucleotide in the sequence is
                    # the key in the frequency_dict
                    if motif[0][i].upper() == key:
                        # Then initialize the position key with value 1
                        # in the value subdictionary
                        value[f"{i}"] = 1
                    # Otherwise, just initialize this subdictionary
                    else:
                        # At count zero
                        value[f"{i}"] = 0
                # Or if the position is already a key in the subdictionary
                elif f"{i}" in value.keys():
                    # Then if the current nucleotide is the same as the
                    # frequency_dict key
                    if motif[0][i].upper() == key:
                        # Then increase the count by 1
                        value[f"{i}"] += 1
    # Save the total number of sequences in this list in the frequency dict
    frequency_dict['total'] = len(motifs_list)
    # And return the frequency dictionary
    return frequency_dict

def calculate_probabilities(frequency_dict):
    """
    Given a frequency dictionary, return a list of lists where
    the zeroeth element is the position of an A or T in the sequence
    and the first element is the probability that the A or T is in
    that position (calculated using Laplace's Definition of Probability)
    """
    assert type(frequency_dict) == dict, "The frequency_dict argument should be a dictionary..."
    needed_keys = ["A", "T", "C", "G", "total"]
    for key in needed_keys:
        assert key in frequency_dict.keys(), f"{key} is not a key in the frequency dictioanry..."
    # Initialize the at_posittions list using list comprehension.
    # The keys of the subdictionary are position strings, and the
    # counts should be initialized to zero
    at_positions = [[key, 0] for key in frequency_dict["A"].keys()]
    # Loop over the position numbers for the sequence, which are the
    # keys of a (any of the) subdictionaries in frequency_dict. This
    # must happen first, as we want the number of As or Ts in the
    # num(th) position of the sequence.  
    for num in frequency_dict["A"].keys():
        # Loop over the keys and values in the frequency_dict
        for key, value in frequency_dict.items():
            # IF the key is total then continue, this is just for probabilities
            if key == "total":
                continue
            # But if the key is A or T,
            elif key == "A" or key == "T":
                # Then loop over the sublists in at_positions
                for sublist in at_positions:
                    # If the position indicator in the sublist is the same
                    # as the num(th) position (the one we're looping over)
                    if sublist[0] == num:
                        # then increase the count for that numth element of the list
                        sublist[1] += value[num]
    # Use list comprehension to calculate the probability of an A or T being
    # in the position.
    at_positions = [[sublist[0], sublist[1]/frequency_dict['total']] for sublist in at_positions]
    # and return these probabilities
    return at_positions

def get_parting_list(at_position_probs, threshold = 0.07):
    """
    Given the AT position probabilities (list of lists) and a threshold
    (default is 0.07), return the positions to partition the motifs on.
    
    The input list to this function is the output from calculate_probabilities
    """
    # Assure that all of the inputs are properly formatted
    assert type(at_position_probs) == list, "The at_position_probs argument should be a list.."
    for pos in at_position_probs:
        assert type(pos) == list, "Each element in the at_position_probs list should also be a list..."
        assert len(pos) == 2, "The subslists should each only be two elements...."
        assert type(pos[0]) == str, "The zeroeth element of each list should be a string..."
        assert type(pos[1]) == float or type(sublist[1]) == int, "The first element of each sublist should be a number..."
    assert type(threshold) == float and threshold >= 0 and threshold <= 0.1, "The threshold should be a float between 0 and 0.1"
    # Initialize the partition_positions list
    partition_positions = []
    # Loop over the positions in the at_position_probs list
    for position in at_position_probs:
        # If the probability of a position is within the thresholded
        # value
        if position[1] <= 0.5+threshold and position[1] >= 0.5-threshold:
            # Then add that position to the list.
            partition_positions.append(int(position[0]))
    # Then return the parition_positions list.
    return partition_positions
            
def get_frequeny_matrix_lines(frequency_dict, tf_name = "TF"):
    """
    Given a frequency dictionary (created by the get_frequency_dict function),
    return a list of lines in JASPAR format to display the frequency
    matrix for the sequence.
    """
    assert type(frequency_dict) == dict, "The frequency dictionary should be a dictionary..."
    assert type(tf_name) == str, "The transcription factor name should be a string.."
    # Initialize the first line of the file, which is the TF name with a carrot
    first_line = [f">{tf_name}\n"]
    # Initialize a list to hold the lines to file
    lines_list = []
    # Loop over the keys and values in the frequency_dict
    for key, value in frequency_dict.items():
        # If the key 'total' shows up, just continue. That stores the total
        # number of sequences
        if key == 'total':
            continue
        # Initialize the line using string formatting. The keys are nucleotides,
        # and the :<3 syntax means "left oriented, allow 3 spaces for the string"
        line = f"{key:<3}["
        # Get the positions list using list comprehension (positions in sequence
        # are keys in each subdictionary in frequency_dict)
        positions = [int(pos) for pos in value.keys()]
        # Sort the positions in ascending order
        positions = sorted(positions)
        # Loop over the positions in the positions list
        for pos in positions:
            # and update the line with the count (subdictionary value at position
            # key), right oriented with 6 characters allowed for the numbers.
            line = f"{line}{value[str(pos)]:>6}"
        # Once the loop ends, add the ending ] character and a newline character
        line = f"{line} ]\n"
        # and add the line to the list
        lines_list.append(line)
    # Once all of the lines are defined, sort the lines based on the first character
    # in each line. 
    lines_list = sorted(lines_list, key = lambda x: x[0])
    # And return the list of the first_line + the lines_list
    return first_line + lines_list


def partition_by_motif_length(motifs_dict):
    """
    Given a motifs dictionary, return a list of the
    motif lists partitioned on length.
    """
    # If only one item is in the motifs dictionary
    if len(motifs_dict.keys()) == 1:
        # Then get the key for that item
        keylist = list(motifs_dict.keys())
        # and return it
        return [motifs_dict[keylist[0]]]
    # Otherwise,
    else:
        # Return them all.
        return [value for key,value in motifs_dict.items()]

def get_partition_sites(sites_file,
                        threshold = 0.07,
                        tf_name = "TF"):
    """
    Given a sites file, a threshold, and a transcription factor name,
    return a list of tuples, where each tuple contains a motif list and
    the equaly likely nucleotides.
    """
    # GEt the directory path to the file
    directory = os.path.split(os.path.realpath(sites_file))[0]
    # Get the motifs from the sites file
    annoted_motifs = parse_sites_file(sites_file)
    # Check the lengths of the motifs
    partitioned_motif_list = check_motif_lengths(annoted_motifs)
    # Partition the motif by length
    partitioned_motif_list = partition_by_motif_length(partitioned_motif_list)
    # Set the motif counter
    motif_counter = 0
    # Initialize the motif list
    motif_list = []
    # loop over the sublists in the paritioned motif list
    for sublist in partitioned_motif_list:
        # Get the frequency dictionary
        frequency = get_frequency_dict(sublist)
        # Calculate the probability of A and T
        prob_list = calculate_probabilities(frequency)
        # get the partition positions based on the probability of A and T
        partition_positions = get_parting_list(prob_list, threshold = threshold)
        # Add the sublist and partition positions to the list
        motif_list.append((sublist, partition_positions))
        # GEt the lines for the jaspar style frequency matrix
        freq_matrix_lines = get_frequeny_matrix_lines(frequency, tf_name)
        # Write the frequency matrix
        write_output_file(freq_matrix_lines, f'{directory}/motif_{motif_counter}.jaspar')
        # Increase the count
        motif_counter += 1
    # Return the motifs list. 
    return motif_list

#
#
#########################################################################
#
#         Functions for making potential binding site lists

def parse_sites_file(a_sites_file,
                     nucleotides_list = nukes):
    """
    Given a .sites file (from JASPAR, download the FASTA version of
    transcription factor binding sites) and a list of nucleotide strings,
    return a list of the TF binding sites.
    
    NOTE: The .sites files are short form FASTA files, where the binding
    sites are all capitalized and the rest of the sequence is not.
    
    The region lines of a JASPAR sites file have the following format:
    >gs_chrNUM:reg_start-reg_end(plus/minus_strand)
    
    an example:
    >dm_chr2L:5980922-5980936(+)
    """
    # Assure that the file name is a string
    assert type(a_sites_file) == str, "The name of the .sites file should be a string"
    # And that it exists
    assert os.path.exists(a_sites_file), f"The file {a_sites_file} does not exist..."
    # Initialize the list to hold the TF binding site strings.
    site_strings = []
    # Open and read the sites file
    with open(a_sites_file, 'r') as f:
        # Initialize a string for a site. This will be updated with
        # each binding site as the file is read.
        site_builder = ""
        # Initialize the go variable. This is used to keep the while loop
        # going, and will be set to False when there are no more lines
        # in the file.
        go = True
        # Start a while loop, using go as the boolean condition for the loop
        while go:
            # Read a line of the file, save it as line
            line = f.readline()
            # Try to get the zeroeth character in the line. Will fail if the
            # line is empty, i.e. the file has no more lines.
            try:
                line[0]
            # IF this does fail, set go to False so the while loop will stop
            except:
                go = False
            # If this does not fail, then check the line type.
            else:
                # If the zeroeth character of the line is a >, you're on a new region.
                # if the site builder has not been initialized, then get the direction
                # of the sequence from the line.
                if line[0] == ">" and site_builder == "":
                    # Then save the strand direction of the upcoming sequence
                    direction = line[-3]
                # Or if the 0th element is a > and the site builder is not an empty string
                elif line[0] == ">" and site_builder != "":
                    # Then save the previuosly built transcription factor and strand direction
                    # as a tuple to the site_strings list
                    site_strings.append((site_builder, direction))
                    # And reset the site_builder variable to the empty string
                    site_builder = ""
                    # and save the strand direction of the next region
                    direction = line[-3]
                # If there is no > character, then we are at a sequence section!
                else:
                    # Make a list of all the nucleotides in the sequence that are also
                    # in the nucleotide list. Nucleotide list is assumed to be capitalized
                    site_list = [letter for letter in line if letter in nucleotides_list]
                    # Update the site_builder variable with these letters, as a string
                    site_builder += str_a_list(site_list)
        # At the end, add the last site string to the list
        site_strings.append((site_builder, direction))
        # Once you have found all of the TF binding sites, the while loop will end and we
        # can close the file.
        f.close()
    # And return the list of TF binding site strings.
    return site_strings     

def minus_to_plus(binding_site_list):
    """
    Given a list of binding sites from a JASPAR sites file,
    return a list of all possible + strand sequences, and all
    possible minus strand sequences. 
    
    NOTE: A sequence on the minus (-) strand will show up as
    the reverse complement on the plus strand. Since the genome
    FASTA files are only the + strand, we need to convert
    all the - strand sequences into + strand sequences. Similarly,
    we can check the minus strands for + strand like sequences by
    finding the reverse complement of the + strand sequences. Therefore,
    all combinations of + and - strand sequences are found.
    
    The revere part is not done in this function, that is saved for
    after partitioning of the sequences.
    """
    assert type(binding_site_list) == list, "The input needs to be a list"
    for site in binding_site_list:
        assert type(site) == tuple, "The input list should only contain strings tuples"
        assert len(site) == 2, "The tuples should only have two elements: (motif_string, +) or (motif_string, -)"
        assert type(site[0]) == str, "The motif should be a string"
        assert "+" in site[1] or "-" in site[1], "The first element of the tuple should be + or -"
    # nucleotide complements in DNA strands
    pairs = {"A" : "T",
             "T" : "A",
             "G" : "C",
             "C" : "G",
             "N" : "N"}
    # Initialize the plus and minus strand lists
    all_possible_plus = []
    all_possible_minus = []
    # Loop over the motifs in the binding_site
    for motif in binding_site_list:
        # Initialize the string_manip variable. The name is a joke
        # from pokemon speedrunning. When they manipulate the computer
        # in some way, they call it a "manip"
        string_manip = ""
        # If the first element of the motif tuple is a +, that
        # means that this motif is on the + strand
        if motif[1] == "+":
            # So add this motif to the + strand list
            all_possible_plus.append(motif[0])
            # and loop over all the letters int he motif
            for letter in motif[0]:
                # and update the string_manip variable with the
                # complement strand
                string_manip = f"{string_manip}{pairs[letter]}"
            # And add the complement to the minus strand list.
            all_possible_minus.append(string_manip)
        # Or if the first element of the motif tuple
        # is a -, then this is a - strand sequence.
        elif motif[1] == "-":
            # So add this motif to the plus strand list
            all_possible_plus.append(motif[0])
            # and loop over all the letters int he motif
            for letter in motif[0]:
                # and update the string_manip variable with the
                # complement strand
                string_manip = f"{string_manip}{pairs[letter]}"
            # And add the complement to the minus strand list.
            all_possible_minus.append(string_manip)
    # At the end, return the lists of all possible plus and minus motifs
    return all_possible_plus, all_possible_minus

def partition_on_unassigned_nucs(motif_list,
                                 partitions_list):
    """
    Given a list of motifs and a list of integers that specify the
    equally likely nucleotide positions in the motifs, return
    a list of lists, where each sublist contains strings split at
    the indicated positions.
    """
    # Assure that the given inputs are the proper types
    assert type(motif_list) == list, "The motifs should be in a list..."
    for motif in motif_list:
        assert type(motif) == str, "Each motif needs to be type str"
    assert type(partitions_list) == list, "The partition positions should be in a list"
    for location in partitions_list:
        assert type(location) == int, "Each partition position should be an integer."
    # Get the size of a motif, minus 1. This is the index of the last element in the
    # motifs list
    size = len(motif_list[0]) - 1
    # If 0 is not in the partitions_list
    if 0 not in partitions_list:
        # Then add zero to the partitions_list
        partitions_list.append(0)
    # If the size of the motif is not in the partitions list,
    if size not in partitions_list:
        # Then add the size to the list as well.
        partitions_list.append(size)
    # And sort the partitions list. This list is used to slice the motifs, so we
    # need to have the beginning, the end, and the partition points in between
    # in order for the following loop to function. Also they need to be in
    # ascending order.
    partitions_list = sorted(partitions_list)
    # Inititalize the parted_motifs list. This will hold all of the lists of 
    # partitioned motifs
    parted_motifs = []
    # Loop over the sequences in the motif list.
    for sequence in motif_list:
        # Initialize the sequence parts list. This will hold the parts of a
        # specific sequence
        seq_parts = []
        # Loop over the number of elements in the partitions list, minus 1.
        for i in range(len(partitions_list)-1):
            # If the ith element of the partitions list is zero
            if partitions_list[i] == 0:
                # Then this we need to slice the string from beginning:first_location
                # and add that substring to the sequence parts list.
                seq_parts.append(sequence[partitions_list[i]:partitions_list[i+1]])
            # Or if the i+1(th) element of the list is the length of the sequence
            elif partitions_list[i+1] == len(sequence)-1:
                # Then we need to slice the string from last_part:end
                seq_parts.append(sequence[partitions_list[i]+1:])
            # If neither the beginning nor the end cases are happening,
            else:
                # Then we are in the middle and we need to slice the string as
                # (location_i)+1 : (location_i+1)-1
                # because location_i needs to be excluded, and so does location_i+1.
                # String slicing tutorial:
                # a = "12345" ->   a[1:3]   ->   "23" (keeps elements 1 and 2, excludes 3)
                seq_parts.append(sequence[partitions_list[i]+1:partitions_list[i+1]])
        # Once the sequence has been partitioned, add the partitioned sequence to the main list
        parted_motifs.append(seq_parts)
    # And at the end, return the parted_motifs list
    return parted_motifs

def bin_the_parts(seqs_split_on_null):
    """
    Given the list of sequences split on the equally likely nucleotides,
    return the parts, binned by their position relative to the original
    sequence.
    """
    # Assure the inputs are in the correct format.
    assert type(seqs_split_on_null) == list, "The input should be a list"
    for item in seqs_split_on_null:
        assert type(item) == list or type(item) == tuple, "The items in the input list should also be lists..."
        for string in item:
            assert type(string) == str, "The input lists should contain lists, which should contain strings......"
    # Use list comprehension to make a list of lists, where each list
    # contains the ith part of split sequences
    bin_list = [[parts[i] for parts in seqs_split_on_null]
                for i in range(len(seqs_split_on_null[0]))]
    # And return the binned list.
    return bin_list

def make_all_combinations(binned_parts,
                          nucleotides_list = nukes):
    """
    Given a list of sequence splits binned and a list of nucleotide strings
    (default is nukes list defined above), generate a combination of a first
    and second sequence part with a new nucleotide between them.
    
    This is a generator function. It does not commit each combination to memory.
    Rather, it generates each combination and returns it, then continues once
    other processes have completed.
    
    This function is kind of slow, but it should create all combinations of
    two lists of strings separated by a new nucleotide.
    """
    # Assure the inputs are all good
    assert type(binned_parts) == list, "The binned_parts argument should be a list..."
    for subbin in binned_parts:
        assert type(subbin) == list, 'The bins in binned_parts should also be lists...'
        for seq in subbin:
            assert type(seq) == str, "Each motif partition should be a str type..."
    # Initialize the current_combo string. This holds the current combination of
    # part_a + N + part_b
    current_combo = ""
    # Initialize the seen list. This will hold all of the combinations that we have
    # already come across, thus avoiding redundancy.
    seen = []
    # Loop over the number of motifs in the zeroeth bin
    for i in range(len(binned_parts[0])):
        # Loop over the number of motifs in the first bin
        for j in range(len(binned_parts[1])):
            # Loop over the nucleotides in the nucleotide list
            for nucleotide in nucleotides_list:
                # If the nucleotide in the list is not N
                if nucleotide.upper() != "N":
                    # Initialize the current combo.
                    current_combo = f"{binned_parts[0][i]}{nucleotide}{binned_parts[1][j]}"
                    # If the current combo is not in the seen list yet,
                    if current_combo not in seen:
                        # Then yield (return and wait).
                        yield current_combo
                        # Once the function resumes, add the current combo to the seen list
                        seen.append(current_combo)
                    # If the current combo is in the seen list
                    else:
                        # Then just continue, we don't need this one.
                        continue

def wrap_all_combinations(binned_parts,
                          nucleotides_list = nukes):
    """
    Given a list of sequence splits binned and a list of nucleotide strings
    (default is nukes list defined above), return a list of all possible
    combinations of the binned parts with all nucleotides.
    
    This function will iterate throught the bins in binned_parts and:
    
    1.) combine bins [0] and [1] using make_all_combinations()
    
    2.) Check if there are other combinations to be made
    
        2a.) If not, then return those combinations.
        2b.) Otherwise, delete bin [0] and make bin[1] equal to the result
             of combining bins [0] and [1], then try again. 
    """
    # Assure the inputs are all good
    assert type(binned_parts) == list, "The binned_parts argument should be a list..."
    for subbin in binned_parts:
        assert type(subbin) == list, 'The bins in binned_parts should also be lists...'
        for seq in subbin:
            assert type(seq) == str, "Each motif partition should be a str type..."
    # Get the number of bins in binned_parts
    looper = len(binned_parts)
    # Loop over the number of bins in binned parts
    for i in range(looper):
        # Use make_all_combinations to generate a list of combinations from the
        # zeroeth and the first bins
        uniq_combos = list( make_all_combinations(binned_parts, nucleotides_list = nukes))
        # If there are only two bins in the binned_parts list
        if binned_parts[2:] == []:
            # Then return the uniq_combos list, you've made all the combinations!!
            return uniq_combos
        # Otherwise,
        else:
            # There are still more things to be combined. First, delete the zeroeth
            # bin as it has already been used.
            del binned_parts[0]
            # Now what was the first bin is the zeroeth bin in binned parts. Since
            # we now need to make all combinations of uniq_combos with binned_parts[1],
            # reassign binned_parts[0] with uniq_combos. 
            binned_parts[0] = uniq_combos
                
def get_reversals(unique_combinations):
    """
    Given a list of uniquely combined motif parts, return the list
    in the reverse order.
    """
    # Assure the input is all good
    assert type(unique_combinations) == list, "The unique_combinations argument must be a list..."
    # Use list comprehension with the reverse_a_str() function to generate alist
    # of all strings in reverse order
    reversals = [reverse_a_str(item) for item in unique_combinations]
    # Return the list of reversed strings.
    return reversals

def combine_reversals_and_clean(unique_combinations,
                                unique_reversals):
    """
    Given two lists, return a combined list of these lists with only unique elements
    """
    # Get the combined list of the two lists
    total = unique_combinations + unique_reversals
    # make the list into a set, and then back into a list.
    total = list(set(total))
    # Return the list
    return total

def format_parted_motifs(motif_list,
                         partition_sites,
                         minus = False):
    """
    Given a motif_list, a list of sites to partition on, and
    whether or not these sequences are for the minus strand
    (default is False), return a list of all combinaitons of
    the motifs.
    """
    assert minus == True or minus == False, "Only boolean values are accepted for named input minus"
    # The assertion statements for each input are included
    # in the functions used in this function.
    #
    # Use partition_on_unassigned_nucs() to get hte parittion list
    parted_list = partition_on_unassigned_nucs(motif_list,
                                               partition_sites)
    # Use bin_the_parts() to bin the parted motifs
    binned_parts = bin_the_parts(parted_list)
    # Use wrap_all_combinations() to get all combinations
    # of the motif parts
    all_motif_combos = wrap_all_combinations(binned_parts)
    # If these motifs are minus strand motifs
    if minus == True:
        # Then use get_reversals to reverse the nucleotide
        # sequences
        all_motif_combos = get_reversals(all_motif_combos)
    return all_motif_combos

#
#
#########################################################################
#
#    Functions for working with peak files

def get_peak_regions(peakfile,
                     extension,
                     delimiter):
    """
    Given the name of a peakfile, the extension for the
    peak file, and a delimiter for the data in the peak
    file, return a dictionary with keys as chromosome regions
    and values as lists of (region_start, region_end) tuples
    """
    # Make sure the file exists and what not
    assert os.path.exists(peakfile), "The given peak file does not exits"
    assert extension in ["narrowPeak", "xls"], "The peak file should be an excel file or a narrowpeak file from MACS3"
    # Use the global peakfile_formats dictionary for formatting
    global peakfile_formats
    # The values in the peakfile_formats dictionary are the columns
    # (math counting) of the desired information. Use list comprehension
    # do get these values (minus 1 for computer counting)
    columns = [int(value)-1 for key,value in peakfile_formats[extension].items()]
    # and sort the columns in ascending order
    columns = sorted(columns)
    # Initialize the chromosome regions dictionary
    chrom_regions = {}
    # Open and read the peak file
    with open(peakfile, 'r') as f:
        # Use list comprehension to get a list of the lines,
        # keeping only the information from the columns
        # of interest.
        try:
            peak_regions = [[line.split(delimiter)[i] for i in columns]for line in f]
        except:
            raise ValueError(f"Invalid delimiter: {delimiter} for this file type")
        # Once the list is attained, close the file
        f.close()
    # Next, loop over the sublists in the peak_regions list
    for sublist in peak_regions:
        # The zeroeth element of the sublist is the chromosome
        # identifier. If this not already in the chrom_regions
        # dictionary
        if sublist[0] not in chrom_regions.keys():
            # Then initialize this region inthe dictionary
            chrom_regions[sublist[0]] = [sublist[1:]]
        # Or if the region is already a key in the dictionary
        elif sublist[0] in chrom_regions.keys():
            # Check to see if the specific region is in the
            # dictionary list already
            if sublist[1:] in chrom_regions[sublist[0]]:
                # If so, then continue
                continue
            # Otherwise
            else:
                # Add this new region to the dictionary
                chrom_regions[sublist[0]].append(sublist[1:])
    # Once the regions have been parsed, return the chrom_regions dictionary
    return chrom_regions


#
#
#########################################################################
#
#    Finding the region sequences in a fasta file

def check_next_region(nuc_count,
                      peak_region_dict,
                      chromosome,
                      region_counter):
    
    """
    Helper function for read_and_parse_fasta()
    
    Given: The current nucleotide count
           A peak_region_dict (dictionary)
           the current chromosome
           the current region within that chromosome
           
    Return: True if the next_region starts in the same region
            thatthe current region ends in
            
            False otherwise
    """
    # No assertions here: These arguments are all passed in
    # from the read_and_parse_fasta() function.
    # The next region index in the list is the current region
    # plus 1
    next_region = region_counter + 1
    # If the next region index is outside of the size of the list
    if next_region >= len(peak_region_dict[chromosome]):
        # Then return False, there is not a region at this index
        return False
    # Or if the next region is the same as the current region
    elif peak_region_dict[chromosome][next_region] == peak_region_dict[chromosome][region_counter]:
        # Then return False
        return False
    # If False hasn't been returned yet, then get the region front
    # and the region end integers
    region_front = int(peak_region_dict[chromosome][next_region][0])
    region_end = int(peak_region_dict[chromosome][next_region][1])
    # If current nucleotide count is between the region parameters
    if nuc_count >= region_front and nuc_count <= region_end:
        # Then return True, we need nucleotides from this region
        return True
    # If not,
    else:
        # Then return False, the function can proceed
        return False

def read_and_parse_fasta(fasta_file,
                         peak_region_dictionary):
    """
    Given the name of a fasta file and a dictionary
    containing peak region information, return
    a dictionary with:
    
    keys: chromosomes/chromosome identifiers
    
    values: lists of tuples:
            (region_start, region_end, fasta_sequence)
    """
    # Assure that the inputs are good
    assert os.path.exists(fasta_file), "The FASTA file given as input does not exist..."
    assert type(peak_region_dictionary) == dict, "The peak_region_dictionary should be of type dict"
    # Initialize the nucleotide count, the region count
    # the fasta_sequence_dictionary, and the iteration count variables.
    nucleotide_count = 0
    region_count = 0
    fasta_sequence_dict = {}
    iter_count = 0
    # Open the fasta file and read
    with open(fasta_file, 'r') as f:
        # Initialize the overlap variable. True tells the function
        # not to move on from the current line and to get more
        # sequences from that line.
        overlap = False
        # Initialize the go variable. This tells the while loop to continue
        go = True
        # Initialize the current chromosome string, the current sequence
        # string and the last region list.
        current_chromosome = ""
        current_sequence = ""
        last_region = []
        # While the go variable is true, continue looping.
        while go:
            # Add one to the iteration count.
            iter_count +=1
            # If the overlap variable is False,
            if overlap == False:
                # Read the line and strip the newline character
                line = f.readline().strip()
            # Next, try to get the zeroeth element of the line
            try:
                line[0]
            # If this fails, then set go to False. The line is empty and the
            # file is therefore empty
            except:
                go = False
            # If the line is not empty, then we need to check the format of the line
            else:
                # If the zeroeth element of the line is a carrot, then this line
                # specifies a region in the FASTA file.
                if line[0] == ">":
                    # Loop over the keys, values inthe peak_region_dictionary
                    for key in peak_region_dictionary.keys():
                        # The keys of the peak_region_dictionary are chromosome
                        # identifiers. Thus, if the key is in the line
                        if key in line:
                            # Then set all of the following variables
                            current_chromosome= key
                            current_sequence = ""
                            region_count = 0
                            nucleotide_count = 0
                            prev_line_len = 0
                            # And initialize the fasta_seque3nce_dictionary with and
                            # empty list at this key
                            fasta_sequence_dict[key] = []
                            # And and break :)
                            break
                # If the line does not have the carrot in it, then it is an actual
                # sequence.
                else:
                    # Take the line and lower all the characters in it.
                    line = line.lower()
                    # If the overlap variable is False (boolean)
                    if overlap == False:
                        # Then increase the nucleotide count by the length of the line
                        # since this is a new line
                        nucleotide_count += len(line)
                    # If the region_count is outside of the length of the list, then
                    # we have exhausted all peak regions in this genomic region.
                    if len(peak_region_dictionary[current_chromosome]) <= region_count:
                        # Therefore, continue until we find a new region
                        continue
                    # Otherwise, we need to analyze the sequences for this genomic region
                    # still.
                    else:
                        # Thus, we should get the region front and region end as integers (for comparisons)
                        region_front = int(peak_region_dictionary[current_chromosome][region_count][0])
                        region_end = int(peak_region_dictionary[current_chromosome][region_count][1])
                        # and update the last region variable.
                        last_region = peak_region_dictionary[current_chromosome][region_count]
                    # Now we have checked the current region we are in, we need to see if the
                    # nucleotide count for this region is within that region or not.
                    if nucleotide_count < region_front:
                        # If the nucleotide count is below the beginning of the region, then continue
                        continue
                    # Or if the nucleotide count is wihtin the current peak region
                    elif nucleotide_count >= region_front and nucleotide_count <= region_end:
                        # Then check to see if the current sequence has been initialized
                        if current_sequence == "":
                            # If it has not, then get the difference between the current
                            # nucleotide count and the beginning of the region. This number
                            # is the distance from the end of the line to the
                            # beginning of the region.
                            dist_from_beg = abs(nucleotide_count - region_front)
                            # Next, check the boundary conditions. If the difference between
                            # the nucleotide count and region front is zero
                            if dist_from_beg == 0:
                                # Then the sequence begins at the begininng of the next line.
                                # so keep the current sequence as zero
                                current_sequence = f""
                            # Or if the difference betweent he nucleotide count adn the beginning
                            # of the region front is the same as the length of the line,
                            elif dist_from_beg == len(line):
                                # Then just initialize the current sequence with the entire line.
                                current_sequence = f"{line}"
                            # In any other case
                            else:
                                # Calculate the distance from the beginning of the line to
                                # the region front
                                dist_to_beg = len(line) - dist_from_beg
                                # And start the current sequence at the nucleotide in that
                                # position.
                                current_sequence = f"{line[dist_to_beg:]}"
                        # If the current sequence is not an empty string,
                        else:
                            # Then simply add the entire line to the current sequence.
                            current_sequence = f"{current_sequence}{line}"
                    # If the nucleotide count is outside of the region_end and the current
                    # sequence is not empty
                    elif nucleotide_count >= region_end and current_sequence != "":
                        # Then calculate the distance from the nucleotide count to the
                        # region end.
                        dist_from_end = abs(nucleotide_count - region_end)
                        # If this distance is zero
                        if dist_from_end == 0:
                            # Then don't add anything. This should have been taken care
                            # of in the previous steps
                            current_sequence= f"{current_sequence}"
                        # In any other case
                        else:
                            # Get the distance from the beginning of the current line to the
                            # end of the region
                            dist_to_beg = len(line) - dist_from_end
                            # And update the current sequence
                            current_sequence= f"{current_sequence}{line[:dist_to_beg]}"
                        # If the current sequence is not empty,
                        if current_sequence != "":
                            # Then update the fasta_sequence_dict with this peak region beginning and end,
                            # as well as the sequence found in this loop
                            fasta_sequence_dict[current_chromosome].append((peak_region_dictionary[current_chromosome][region_count][0],
                                                                            peak_region_dictionary[current_chromosome][region_count][1],
                                                                            current_sequence))
                            # Then reset the current sequence
                            current_sequence = ""
                            # and increase the region count by 1
                            region_count += 1
                            # Check for overlap between this region and the next region
                            overlap = check_next_region(nucleotide_count,
                                                        peak_region_dictionary,
                                                        current_chromosome,
                                                        region_count)
                            
    # At the end of this huge loop, return the sequence dictionary
    return fasta_sequence_dict


def check_parsed_fastas(fasta_seq_dict):
    """
    Given a fasta_sequence_dictionary, check to make sure the sequence
    lengths match the theoretical length of the region and toss
    the sequences that do not match.
    """
    # Assure the inputs are formatted properly
    assert type(fasta_seq_dict) == dict, "The fasta sequence dictionary should have type dict..."
    for key, value in fasta_seq_dict.items():
        assert type(value) == list, "The dictionary should have list values..."
        for item in value:
            assert type(item) == tuple, "The lists should all have tuples in them..."
            assert len(item) == 3, "The tuples should all have 3 elements..."
    
    # Initialize the variables for total number of sequences found,
    # total number of sequences kept and total tossed
    total = 0
    total_kept = 0
    total_tossed = 0
    # Initialize a list of keys (chromosome regions) that are
    # tossed because they are empty
    keys_to_toss = []
    # Loop over the keys and values in the fasta_seq_dict
    for key, value in fasta_seq_dict.items():
        # Initialize a list of indices to remove, if they fail
        # the test
        remove_indices = []
        # Loop over each sublist in the current dict value
        for sublist in value:
            # If the length of the sequence and the theoretical distance
            # are different, then tell the user
            if int(sublist[1]) - int(sublist[0]) != len(sublist[2]):
                print(f"Something exploded while extracting the following FASTA sequence:\n")
                print(f"Chromosome: {key}\t Sequence Start: {sublist[0]}\t Sequence End: {sublist[1]}")
                print(f"Extracted Sequence Length: {len(sublist[2])}\t Theoretical Length: {int(sublist[1]) - int(sublist[0])}\n")
                print(f"This sequence will be omitted for the remainder of the program. \n")
                # Get the index for these values
                remove_indices.append(value.index(sublist))
                # And increase the tossed count by one
                total_tossed += 1
            # Otherwise
            else:
                # Increase the kept count by 1
                total_kept += 1
            # Regardless, increase the total count by one
            total += 1
        # Once this has ceased, sort the removal indices in
        # reverse order
        remove_indices= sorted(remove_indices, reverse = True)
        # Loop over the indices
        for index in remove_indices:
            # And remove them in reverse order
            del value[index]
            if fasta_seq_dict[key] == []:
                # If the value is now empty, add the key to
                # the keys_to_toss list
                keys_to_toss.append(key)
    # Loop over the keys and values in the keys_to_toss
    for key in keys_to_toss:
        # and delete those keys
        del fasta_seq_dict[key]
    # Then tell the users the results
    print(f"FASTA sequences checked for: Sequence Length matching the Theoretical Length.")
    print(f"Total sequences found:   {total}")
    print(f"Total sequences kept:    {total_kept}")
    print(f"Total sequences tossed:  {total_tossed}\n")
    if keys_to_toss != []:
        print(f"The following keys were removed because tossing sequences resulted in empty lists:")
        
        for key in keys_to_toss:
            print(f"{key}")
    print("Proceeding with the remainder of the program :)")

#
#
#########################################################################
#
#

def true_false_fastas_with_binding_seqs(possible_binding_motifs,
                                        fasta_seq_dict):
    """
    Given a list of possible binding motifs and a fasta sequence
    dictionary, return two dictionaries:
    
    The true_dict holds all sequences for which there was a match
    between a possible binding motif and a fasta sequence dictionary
    
    The false_dict holds all sequences for which no match was found.
    """
    # Assure that all of the inputs are formatted properly
    assert type(possible_binding_motifs) == list, "Binding motifs should be in a list"
    for item in possible_binding_motifs:
        assert type(item) == str, "Binding motifs should be strings... How did you break this?"
    assert type(fasta_seq_dict) == dict, "The fasta sequences should be in a dictionary"
    # Initialize the true_dict and false_dict
    true_dict = {}
    false_dict = {}
    # Loop over the keys and values in the fasta_seq_dict
    for key, value in fasta_seq_dict.items():
        # Loop over the number of items in the value list
        for i in range(len(value)):
            # Initialize the updated_sequence and location variables
            updated_sequence = ""
            location = ""
            # Loop over the motifs in the possible_binding_motifs list
            for motif in possible_binding_motifs:
                # IF there is no updated sequence and the motif is in the sequence
                if updated_sequence == "" and motif.lower() in value[i][2]:
                    # Then update the sequence and save it as the updated sequence
                    updated_sequence = value[i][2].replace(motif.lower(), motif)
                    # If the location is also empty
                    if location == "":
                        # Then update the location string with the motif and where the motif
                        # is found in the sequence
                        location = f"{motif}:({updated_sequence.index(motif) + 1 + int(value[i][0])})"
                    # Otherwise, update the location
                    else:
                        # With the new information
                        location = f"{location},{motif}:({updated_sequence.index(motif) + 1 + int(value[i][0])})"
                # Or if the updated sequence has been identified and another motif is in that sequence
                elif updated_sequence != "" and motif.lower() in updated_sequence:
                    # Then again update the sequence
                    updated_sequence = updated_sequence.replace(motif.lower(), motif)
                    # And if the location is (for some reason) empty, update it
                    if location == "":
                        location = f"{motif}:({updated_sequence.index(motif) + 1 + int(value[i][0])})"
                    # Or if the location is not empty, add this new location to it
                    else:
                        location = f"{location},{motif}:({updated_sequence.index(motif) + 1 + int(value[i][0])})"
                # If the motif is not in the sequence, then continue
                else:
                    continue
            # If the location and the updated_sequence variables are empty
            if location =="" and updated_sequence == "":
                # Then this goes in the false dictionary. Check to see if the key
                # already exists in the false dictionary
                if key not in false_dict.keys():
                    # If not, then initialize the value
                    false_dict[key] = [(value[i][0], value[i][1], value[i][2], None, False)]
                # If so, then simply add the value
                else:
                    false_dict[key].append((value[i][0], value[i][1], value[i][2], None, False))
            # And if th elovation and updated sequence are not empty
            else:
                # Then this goes in the true dictionary. Again, check to see if the key
                # is already in the dictionary
                if key not in true_dict.keys():
                    # If not, initialize this list
                    true_dict[key] = [(value[i][0], value[i][1], updated_sequence, location, True)]
                # And if so, add this sequence to the list.
                else:
                    true_dict[key].append((value[i][0], value[i][1], updated_sequence, location, True))
    # In the end, return the true and false dictionaries.
    return true_dict, false_dict

#
#
#########################################################################
#
#        Functions for generating lines for file writing

def fasta_dict_to_lines(fasta_seq_dict,
                        tf_name=None):
    """
    Given a fasta sequence dictionary and a transcription factor
    name (default None), return a list of lines formatted
    for FASTA files
    """
    # If the fasta_seq_dict is empty, then return an empty list
    if fasta_seq_dict == {}:
        return []
    # Otherwise, initialize the lines list
    lines = []
    # Loop over the keys and values in the fasta sequence dict
    for key, value in fasta_seq_dict.items():
        # loop over the number of items in the value list
        for i in range(len(value)):
            # If the sequence in the ith element of value is
            # empty, then continue
            if value[i][2] == "":
                continue
            # Otherwise
            else:
                # Add a newline character every 80 characters
                sequence = add_newline_every_x_chars(value[i][2])
                # If there was no TF name given, then just make the
                # line have the chrom, the region start and end
                if tf_name == None:
                    newline = f">{key}\tregion:{value[i][0]}-{value[i][1]}\n{sequence}\n"
                # If a TF name was given, then use that information
                # in the line generation
                else:
                    newline = f">{key}\tregion:{value[i][0]}-{value[i][1]}\tpossible {tf_name} binding elements:{value[i][3]}\n{sequence}\n"
                # and add this string to the lines list
                lines.append(newline)
    # At the end, return the lines list.
    return lines

def combine_true_false_seq_lines(true_fasta_dict,
                                 false_fasta_dict,
                                 tf_name="TF"):
    """
    Given a true and false fasta dictionary pair and a transcription
    factor name, return a list of lines to write with the
    true lines first
    """
    assert type(true_fasta_dict) == dict, "true_fasta_dict should be of type 'dict'"
    assert type(false_fasta_dict) == dict, "false_fasta_dict should be of type 'dict'"
    # Use the fasta_dict_to_lines() to get true lines
    true_lines = fasta_dict_to_lines(true_fasta_dict, tf_name)
    # Use the fasta_dict_to_lines() to get false lines
    false_lines = fasta_dict_to_lines(false_fasta_dict, tf_name)
    # If both lists are empty, raise a type error and exit
    if true_lines == [] and false_lines == []:
        raise TypeError(f"Something exploded, and no sequences are available to write to a file..")
    # Otherwise, combine the lists and return that list
    else:
        lines_to_write = true_lines + false_lines
        return lines_to_write

#
#
#########################################################################
#
#

def write_output_file(lines_to_write,
                      output_filename = "peak_sequences.fasta"):
    """
    Given a list of lines to write and an output file name,
    write a file!
    """
    # Open the file with the writing tag
    with open(output_filename, 'w') as f:
        # Use the writelines method to write the lines
        f.writelines(lines_to_write)
        # close the file
        f.close()

#
#
#########################################################################
#
#        Getting & Writing a Fasta file with binding motifs
    
def fasta_with_binding_motifs(sites_file, 
                              peak_file,
                              genome_fasta,
                              output_filename = "peak_sequences.fasta",
                              use_partition = False,
                              tf_name = "TF",
                              delimiter = '\t'):
    print('Beginning:\n')
    print(f'Getting the annotated motifs from the file {sites_file}...\n')
    print(f"Checking the annotated motifs for partition sites...\n")
    annoted_motifs = get_partition_sites(sites_file, tf_name = tf_name)
    assert len(annoted_motifs) == 1, "Motifs of different lengths are present in your sites file... Please only have one motif length in a file..."
    print(f"Getting all Plus and Minus strand sequences\n")
    all_plus, all_minus = minus_to_plus(annoted_motifs[0][0])
    if annoted_motifs[0][1] != [] and use_partition == True:
        print(f"An equally likely nucleotide position was found! Getting all possible sequences...\n")
        all_plus = format_parted_motifs(all_plus, annoted_motifs[0][1], minus = False)
        all_minus = format_parted_motifs(all_minus, annoted_motifs[0][1], minus = True)
        all_motifs = combine_reversals_and_clean(all_plus,
                                                 all_minus)
        del all_plus
        del all_minus
        del annoted_motifs
    else:
        print(f"No equally likely nucleotides were found. Proceeding using original sequences...\n")
        all_minus = get_reversals(all_minus)
        all_motifs = combine_reversals_and_clean(all_plus,
                                                 all_minus)
        del all_plus
        del all_minus
        del annoted_motifs
    peak_extension = peak_file.split('.')[-1]
    print(f'Gathering peak regions from {peak_file}...\n')
    peak_regions_dict = get_peak_regions(peak_file, peak_extension, delimiter)
    print(f"Finding the genomic sequences corresponding to peak regions...\n")
    peak_sequence_dict = read_and_parse_fasta(genome_fasta, peak_regions_dict)
    del peak_regions_dict
    print(f"Checking the sequences for consistency...\n")
    check_parsed_fastas(peak_sequence_dict)
    print(f"\nMapping motifs to the sequences found...\n")
    true_sites_dict, false_sites_dict = true_false_fastas_with_binding_seqs(all_motifs,
                                                                            peak_sequence_dict)
    del all_motifs
    del peak_sequence_dict
    print(f"Creating the lines for the FASTA file...\n")
    lines = combine_true_false_seq_lines(true_sites_dict,
                                         false_sites_dict,
                                         tf_name=tf_name)
    del true_sites_dict
    del false_sites_dict
    print(f"Writing the FASTA file...\n")
    write_output_file(lines, output_filename=f"{os.path.split(os.path.realpath(peak_file))[0]}/{output_filename}")
    del lines
    print(f"Done! The file can be found in {os.path.dirname(peak_file)}...")



#
#
#########################################################################
#
#
    
def fasta_without_binding_motifs(peak_file,
                                 genome_file,
                                 output_filename = "peak_sequences.fasta",
                                 delimiter = '\t'):
    """
    Given a peak file, a genome_file (fasta format), an output filename, and delimiter,
    find the peak region sequences and write them to a fasta file.
    """
    peak_extension = peak_file.split('.')[-1]
    peak_regions_dict = get_peak_regions(peak_file, peak_extension, delimiter)
    peak_sequence_dict = read_and_parse_fasta(genome_fasta, peak_regions_dict)
    del peak_regions_dict
    check_parsed_fastas(peak_sequence_dict)
    lines = fasta_dict_to_lines(peak_sequence_dict, tf_name = None)
    del peak_sequence_dict
    write_output_file(lines, output_filename = output_filename)
    del lines
    print('done')



#
#
#########################################################################
#
#