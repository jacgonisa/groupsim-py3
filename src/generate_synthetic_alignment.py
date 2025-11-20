#!/usr/bin/env python
"""
generate_synthetic_alignment.py - Creates a simple FASTA alignment for testing GroupSim.
"""

import random

def generate_synthetic_alignment(num_seqs=10, length=30, output_prefix="synthetic_test"):
    
    # 5 sequences in Group 1, 5 sequences in Group 2
    num_group1 = num_seqs // 2
    num_group2 = num_seqs - num_group1
    
    # Define residues for control columns
    CONSERVED_RES = 'A'
    SDP1_RES = 'K'
    SDP2_RES = 'D'
    
    # General amino acid set for random columns
    AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
    
    # Sequence list, headers list
    sequences = []
    headers = []
    
    # --- Generate Sequences ---
    
    for i in range(num_seqs):
        
        is_group1 = (i < num_group1)
        group_id = 1 if is_group1 else 2
        header = f"Seq{i+1}|Group{group_id}"
        headers.append(header)
        
        sequence = []
        
        for pos in range(length):
            if pos == 5:
                # Position 5: Perfectly conserved (Control) -> Score near zero
                sequence.append(CONSERVED_RES) 
            elif pos == 15:
                # Position 15: Perfect SDP (SDP) -> Max positive score
                sequence.append(SDP1_RES if is_group1 else SDP2_RES)
            elif pos == 25:
                # Position 25: Variable within groups, different between (Moderate SDP)
                group_res = ('R', 'K', 'H') if is_group1 else ('E', 'D', 'Q')
                sequence.append(random.choice(group_res))
            else:
                # Other positions: Random noise
                sequence.append(random.choice(AMINO_ACIDS))
                
        sequences.append("".join(sequence))

    # --- Write Files ---
    
    # 1. FASTA Alignment
    fasta_file = f"{output_prefix}.fasta"
    with open(fasta_file, 'w') as f:
        for header, seq in zip(headers, sequences):
            f.write(f">{header}\n{seq}\n")
            
    # 2. Manual Group File
    group_file = f"{output_prefix}_groups.txt"
    with open(group_file, 'w') as f:
        f.write("# Manual groups for synthetic alignment\n")
        f.write("Group1: " + ", ".join(h.split('|')[0] for h in headers[:num_group1]) + "\n")
        f.write("Group2: " + ", ".join(h.split('|')[0] for h in headers[num_group1:]) + "\n")

    print(f"Generated files:\n- {fasta_file}\n- {group_file}")
    print("\nTo run GroupSim, use: ")
    print(f"python auto_groupsim_final_v3.py -k {group_file} {fasta_file}")

if __name__ == "__main__":
    generate_synthetic_alignment()
