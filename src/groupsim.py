#!/usr/bin/env python
"""
groupsim.py - Automated/Manual GroupSim with Full Visualization and Enhanced Plotting

Features:
1. Automatic Clustering by Identity Cutoff (-t).
2. Automatic Clustering by Target Number of Groups (-k N).
3. Manual Group Definition via file (-k filename) - Unlisted sequences become 'Other'.
4. Visualization: Clustered Heatmap, Dendrogram, and GroupSim Score Plot.
5. Score Plot Annotation: Includes position and residue consensus (e.g., G1:A/T | Other:S).
6. Enhanced Score Plot Legend: Uses a gradient color bar for Z-scores.

Requires: biopython, numpy, pandas, scipy, matplotlib, seaborn
"""

import sys
import math
import getopt
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import AlignIO
from scipy.spatial.distance import squareform
from scipy.cluster import hierarchy
from scipy.stats import zscore


# --- Global Definitions ---
amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
	       'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
amino_acids_nogap = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
		     'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

# Identity Matrix for default scoring
ident = []
for i in range(21):
    t = []
    for j in range(21):
	    t.append(1.0 if i == j else 0.0)
    ident.append(t)


def usage():
    """Prints the usage message for the script."""
    usage_string = f"""
USAGE:
auto_groupsim_final_v3.py [options] alignfile

    -You must specify EITHER -t or -k.

OPTIONS for Group Inference:
    -t [real in [0, 100]] (AUTOMATIC MODE)
     Identity Cutoff. Groups sequences closer than this percentage (e.g., -t 75.0).

    -k [int >= 2] (AUTOMATIC MODE)
     Target Number of Groups (K). Finds the necessary cutoff to create K clusters (e.g., -k 3).
     
    -k [filename] (MANUAL MODE)
     File containing manual group definitions. Format: 'Group1:SeqA,SeqB' (one line per group/sequences)
     Sequences not listed are placed in the 'Other' group.

    -o [filename_prefix]
     Name prefix for output files. Default=alignment filename prefix.

OPTIONS for GroupSim Scoring:
    -c [real in [0, 1)] (Column Gap Cutoff. Default=.1)
    -g [real in [0, 1)] (Group Gap Cutoff. Default=.3)
    -l [real in [0,1]] (Lambda for window heuristic. Default=.7)
    -m [filename] (Similarity matrix file. Default=identity matrix)
    -n (Map raw scores to [0,1].)
    -w [int] (Conservation window size. Default=3)
    -h (Help. Print this message.)
"""
    print(usage_string, file=sys.stderr)

# --- Core Utility Functions ---

def read_scoring_matrix(sm_file):
    """Read a scoring matrix from a file and return it."""
    if not sm_file: return ident, 'identity'
    try:
        with open(sm_file, 'r') as f:
            lines = [line.strip().split() for line in f if not line.startswith('#') and line.strip()]
        
        if len(lines) >= 20 and len(lines[0]) < 20: # Handle triangular matrix format
            for i in range(20):
                for j in range(i + 1, 20):
                    if j < len(lines) and i < len(lines[j]): lines[i].append(lines[j][i])
        
        matrix = [[float(val) for val in row] for row in lines if len(row) >= 20]
        return matrix, sm_file
    except Exception as e:
        print(f"Error loading matrix {sm_file}: {e}. Using identity matrix...", file=sys.stderr)
        return ident, 'identity'

def read_alignment(filename):
    """Reads alignment in FASTA or CLUSTAL format using BioPython."""
    try:
        alignment = AlignIO.read(filename, "clustal")
    except:
        try:
            alignment = AlignIO.read(filename, "fasta")
        except Exception as e:
            sys.exit(f"Problem reading or parsing alignment file {filename}: {e}. Exiting...")
            
    names = [record.id for record in alignment]
    sequences = [str(record.seq).upper().replace('B', 'D').replace('Z', 'Q').replace('X', '-') for record in alignment]
        
    if not sequences or not all(s for s in sequences): sys.exit("Alignment is empty or sequences could not be parsed.")
    return names, sequences

def get_column(col_num, alignment):
    """Return the col_num column of alignment as a list."""
    return [seq[col_num] for seq in alignment]

def format_column(col_num, alignment, gids_to_seqs):
    """Return a string with residues for each group in col_num column."""
    col = get_column(col_num, alignment)
    groups = []
    
    for gid in sorted(gids_to_seqs.keys()):
        groups.append("".join(col[seq_index] for seq_index in gids_to_seqs[gid]))

    return " | ".join(groups)

def gap_percentage(col):
    """Return the percentage of gaps in col."""
    return col.count('-') / len(col) if col else 0.0

def norm_01(score_list):
    """Map the values in score_list to the range [0,1]. Ignore None."""
    valid_scores = [s for s in score_list if s is not None]
    if not valid_scores: return score_list

    mini, maxi = min(valid_scores), max(valid_scores)
    range_val = float(maxi - mini)
    
    new_sl = []
    for s in score_list:
        if s is None:
            new_sl.append(None)
        else:
            new_sl.append(1.0 if range_val == 0 else (s - mini) / range_val)

    return new_sl

def freq_count(col, symbols, pc=0.000001):
    freq_counts = {sym: pc for sym in symbols}
    for sym in col:
        if sym in freq_counts: freq_counts[sym] += 1
    total_counts = float(sum(freq_counts.values()))
    if total_counts == 0: return [0.0] * len(symbols)
    return [freq_counts[sym] / total_counts for sym in symbols]

def js_divergence(col1, col2=None):
    fc1 = freq_count(col1, amino_acids_nogap, 0.001)
    # Background distribution (BLOSUM62 approximation) if col2 is None
    fc2 = [0.074, 0.052, 0.045, 0.054, 0.025, 0.034, 0.054, 0.074, 0.026, 0.068, 0.099, 0.058, 0.025, 0.047, 0.039, 0.057, 0.051, 0.013, 0.032, 0.073] if not col2 else freq_count(col2, amino_acids_nogap, 0.001)

    if len(fc1) != len(fc2): return -1
    r = [0.5 * fc1[i] + 0.5 * fc2[i] for i in range(len(fc1))]
    d = 0.0
    for i in range(len(fc1)):
        term = 0.0
        if fc1[i] > 0.0 and r[i] > 0.0: term += fc1[i] * math.log(fc1[i] / r[i], 2)
        if fc2[i] > 0.0 and r[i] > 0.0: term += fc2[i] * math.log(fc2[i] / r[i], 2)
        d += term

    return (1.0 - gap_percentage(col1)) * (d / 2)

def cons_window_score(sdp_scores, alignment, window_len, lam):
    win_scores = sdp_scores[:]
    columns = [get_column(i, alignment) for i in range(len(alignment[0]))]
    cons_scores = [js_divergence(col) if sdp_scores[i] is not None else None for i, col in enumerate(columns)]
    wincon_scores = []

    for i, sdp_score in enumerate(sdp_scores):
        if sdp_score is None:
            wincon_scores.append(None)
            continue
        sum_con, num_terms = 0.0, 0.0
        start, end = max(i - window_len, 0), min(i + window_len + 1, len(sdp_scores))
        
        for j in range(start, end):
            if i != j and cons_scores[j] is not None and cons_scores[j] >= 0:
                num_terms += 1.0
                sum_con += cons_scores[j]
                
        wincon_scores.append(sum_con / num_terms if num_terms != 0 else 0.0)

    wincon_scores = norm_01(wincon_scores)
    for i, sdp_score in enumerate(sdp_scores):
        if sdp_score is not None and wincon_scores[i] is not None:
            win_scores[i] = (1.0 - lam) * wincon_scores[i] + lam * sdp_score
    return win_scores

def matrix_sum_pairs(groups, matrix=ident):
    sum_within_groups, num_within_groups = 0.0, 0
    
    for group in groups:
        local_sum, num_terms = 0.0, 0
        for i in range(len(group)):
            for j in range(i):
                try:
                    i1, i2 = amino_acids.index(group[i]), amino_acids.index(group[j])
                    if i1 < 20 and i2 < 20: 
                         local_sum += matrix[i1][i2]
                         num_terms += 1
                except ValueError: pass

        if num_terms > 0: 
            sum_within_groups += (local_sum / num_terms)
            num_within_groups += 1

    mean_within = sum_within_groups / num_within_groups if num_within_groups > 0 else 0.0

    sum_between_groups, num_between_pairs = 0.0, 0
    for g1_ind in range(len(groups)):
        for g2_ind in range(g1_ind):
            pair_sum, num_terms = 0.0, 0
            for aa1 in groups[g1_ind]:
                for aa2 in groups[g2_ind]:
                    try:
                        i1, i2 = amino_acids.index(aa1), amino_acids.index(aa2)
                        if i1 < 20 and i2 < 20:
                            pair_sum += matrix[i1][i2]
                            num_terms += 1
                    except ValueError: pass

            if num_terms > 0:
                sum_between_groups += (pair_sum / num_terms)
                num_between_pairs += 1

    mean_between = sum_between_groups / num_between_pairs if num_between_pairs > 0 else 0.0
    group_residues_flat = "".join("".join(g) for g in groups)
    conservation_weight = (1.0 - gap_percentage(group_residues_flat))
    
    return conservation_weight * (mean_within - mean_between)

# --- Group Inference Functions ---

def calculate_similarity_matrix(names, alignment):
    """Calculates the pairwise percentage identity matrix and returns it as a DataFrame."""
    n = len(names)
    sequences = alignment 
    similarity_matrix = np.zeros((n, n))
    
    for i in range(n):
        seq_i = sequences[i]
        for j in range(i, n):
            seq_j = sequences[j]
            matches = sum(a == b for a, b in zip(seq_i, seq_j) if a != "-" and b != "-")
            
            len_i_nogap = len(seq_i.replace("-", ""))
            len_j_nogap = len(seq_j.replace("-", ""))
            length = min(len_i_nogap, len_j_nogap)
            
            similarity = (matches / length) * 100.0 if length > 0 else 0.0
            similarity_matrix[i, j] = similarity
            similarity_matrix[j, i] = similarity
            
    return pd.DataFrame(similarity_matrix, index=names, columns=names)

def parse_manual_groups_from_file(filename, names):
    """
    Parses a group definition file and assigns unlisted sequences to 'Other'.
    File format: 'Group1:SeqA,SeqB' (one line per group/sequences)
    """
    name_to_index = {name: i for i, name in enumerate(names)}
    assigned_indices = set()
    gids_to_seqs = {}
    
    print(f"Reading manual groups from file: {filename}", file=sys.stderr)
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'): continue
                
                if ':' not in line: continue
                    
                gid, seq_names_str = line.split(':', 1)
                gid = gid.strip()
                # Split sequence names by comma and/or space
                seq_names = [s.strip() for s in seq_names_str.replace(',', ' ').split() if s.strip()]
                
                if not gid: continue
                if gid not in gids_to_seqs: gids_to_seqs[gid] = []
                
                for seq_name in seq_names:
                    if seq_name in name_to_index:
                        seq_index = name_to_index[seq_name]
                        if seq_index not in assigned_indices:
                            gids_to_seqs[gid].append(seq_index)
                            assigned_indices.add(seq_index)
                    else:
                        print(f"WARNING: Sequence '{seq_name}' not found in alignment. Skipping.", file=sys.stderr)
    except IOError:
        sys.exit(f"ERROR: Could not open group file {filename}.")

    # Assign all unlisted sequences to the default 'Other' group
    other_group_id = 'Other'
    gids_to_seqs[other_group_id] = [i for i in range(len(names)) if i not in assigned_indices]

    if len(gids_to_seqs) < 2:
        sys.exit("ERROR: Manual grouping resulted in less than two groups (even with the default 'Other' group).")
        
    print(f"Created {len(gids_to_seqs)} groups: {', '.join(gids_to_seqs.keys())}.", file=sys.stderr)
    return gids_to_seqs, 100.0 # 100% identity is used as a placeholder for plotting

# --- Visualization Functions ---

def plot_clustered_heatmap(sim_df, linkage_matrix, output_name):
    """Generates and saves a clustered heatmap of the similarity matrix."""
    fig_name = f"{output_name}_heatmap.png"
    g = sns.clustermap(
        sim_df, row_linkage=linkage_matrix, col_linkage=linkage_matrix,
        cmap="viridis", figsize=(12, 12), cbar_kws={'label': 'Percent Identity'},
        linewidths=0.5, linecolor="black", dendrogram_ratio=(.1, .2) 
    )
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=6)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=6)
    g.fig.suptitle("Hierarchical Clustering of Sequence Similarities", y=1.02)
    plt.savefig(fig_name, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Wrote clustered heatmap to {fig_name}.", file=sys.stderr)

def plot_dendrogram(linkage_matrix, names, identity_cutoff, output_name):
    """Generates and saves the dendrogram with the group cutoff line."""
    fig_name = f"{output_name}_dendrogram.png"
    plt.figure(figsize=(10, max(6, len(names) * 0.3))) 
    distance_cutoff = 1.0 - (identity_cutoff / 100.0)
    
    hierarchy.dendrogram(
        linkage_matrix, orientation='right', labels=names, distance_sort='descending',
        show_leaf_counts=True, color_threshold=distance_cutoff, leaf_font_size=8,
    )
    
    plt.axvline(x=distance_cutoff, color='r', linestyle='--', 
                label=f'{identity_cutoff:.1f}% Identity Cutoff')
    
    def dist_to_id(x): return 100.0 * (1.0 - x)
        
    sec_ax = plt.gca().secondary_xaxis('top', functions=(dist_to_id, dist_to_id))
    sec_ax.set_xlabel('Percent Identity')
    
    plt.xlabel('Distance (1 - Identity/100)')
    plt.ylabel('Sequence ID')
    plt.title(f'Hierarchical Clustering Dendrogram (Cutoff: {identity_cutoff:.1f}%)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(fig_name, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Wrote dendrogram to {fig_name}.", file=sys.stderr)

def get_residue_consensus(col_residues):
    """Generates a string representing the consensus residues per group (e.g., G1:A/T | G2:S)."""
    
    consensus_parts = []
    for gid in sorted(col_residues.keys()):
        residues = [r for r in col_residues[gid] if r != '-']
        
        if not residues:
            consensus_parts.append(f"{gid}:-")
        else:
            unique_residues = sorted(list(set(residues)))
            # Concisely represent the consensus: e.g., A/T/S
            consensus_parts.append(f"{gid}:{'/'.join(unique_residues)}")
            
    return " | ".join(consensus_parts)


def plot_groupsim_manhattan(scores, alignment, gids_to_seqs, output_name, threshold=2.0):
    """
    Generates and saves a GroupSim score plot showing SDP scores across alignment positions,
    annotating sites with the residue change, and using a gradient Z-score legend.
    """
    
    fig_name = f"{output_name}_manhattan_plot.png"
    valid_scores = [(i, score) for i, score in enumerate(scores) if score is not None]
    if not valid_scores:
        print("No valid scores to plot.", file=sys.stderr); return
        
    positions, score_values = zip(*valid_scores)
    score_array = np.array(score_values)
    mean_score, std_score = np.mean(score_array), np.std(score_array)
    
    z_scores = (score_array - mean_score) / std_score if std_score != 0 else np.zeros_like(score_array)

    plot_df = pd.DataFrame({
        'Alignment position': positions, 'Groupsim score': score_values, 'z score': z_scores
    })

    plt.figure(figsize=(12, 6))
    
    # Scatter Plot (colored by Z-score with a color bar)
    ax = sns.scatterplot(
        x='Alignment position', y='Groupsim score', hue='z score', data=plot_df, 
        palette='coolwarm', hue_norm=(0, max(10, z_scores.max() if z_scores.size > 0 else 0)), 
        s=50, legend="full",  # "full" to show the color bar
        edgecolor='black', linewidth=0.5
    )
    
    # Add a title to the color bar
    norm = plt.Normalize(vmin=0, vmax=max(10, z_scores.max() if z_scores.size > 0 else 0))
    sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=norm)
    sm.set_array([]) # dummy array for the color bar
    ax.figure.colorbar(sm, ax=ax, label="Z-score")

    # Line Plot (Smoothed trend/mean line) - Make it appear on top
    sns.lineplot(
        x='Alignment position', y='Groupsim score', data=plot_df, 
        color='black', linewidth=2, alpha=0.8, estimator='mean', ci=None, zorder=10 # Higher zorder to be on top
    )
    
    plt.title("GroupSim Scores Across Alignment Positions", fontsize=14, fontweight='bold')
    plt.xlabel("Alignment position", fontsize=12)
    plt.ylabel("Groupsim score", fontsize=12)
    
    # Labeling significant scores (Z > threshold)
    significant_sites = plot_df[plot_df['z score'] > threshold]
    
    if not significant_sites.empty:
        print(f"Outlier sites (Z > {threshold}):", file=sys.stderr)
        
        for index, row in significant_sites.iterrows():
            pos = int(row['Alignment position'])
            score = row['Groupsim score']
            
            # Get residues for this column, grouped by GID
            col_residues = {}
            for gid in gids_to_seqs:
                col_residues[gid] = [alignment[seq_idx][pos] for seq_idx in gids_to_seqs[gid]]
            
            residue_change = get_residue_consensus(col_residues)
            
            # Annotation Text: Use only the residue change, pos is on X-axis
            ann_text = residue_change
            
            plt.annotate(
                ann_text,
                (pos, score), textcoords="offset points", xytext=(0, 10), ha='center',
                fontsize=8, arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0.1", color='gray')
            )
            print(f"  Pos {pos}: Score={score:.3f}, Residues: {residue_change}", file=sys.stderr)

    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(fig_name, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Wrote GroupSim score plot to {fig_name}.", file=sys.stderr)

def write_detailed_output(score_file, names, alignment, scores, gids_to_seqs, group_ids, final_identity_cutoff):
    """Writes a detailed, column-by-column report."""
    
    index_to_gid = {}
    for gid in gids_to_seqs:
        for index in gids_to_seqs[gid]: index_to_gid[index] = gid
    
    score_file.write("\n## Detailed Specificity Score and Alignment\n")
    score_file.write(f"# Grouping Cutoff: {final_identity_cutoff:.1f}% identity\n")

    num_cols = len(alignment[0]) if alignment else 0
    max_name_len = max(len(name) + len(index_to_gid.get(i, 'N/A')) + 1 for i, name in enumerate(names)) if names else 15
    max_name_len = max(max_name_len, 20) 

    score_line = f"{'Score':<{max_name_len}}"
    index_line = f"{'Index':<{max_name_len}}"
    
    for i in range(num_cols):
        score_str = f'{scores[i]:.3f}' if scores[i] is not None else 'None '
        score_line += f'{score_str:<7}'
        index_line += f'{i:<7}'
    
    score_file.write(score_line.rstrip() + '\n')
    score_file.write(index_line.rstrip() + '\n')
    score_file.write('-' * (max_name_len + num_cols * 7) + '\n')
    
    for i, name in enumerate(names):
        gid = index_to_gid.get(i, 'N/A')
        label = f"{name}|{gid}"
        residues = "".join(alignment[i])
        
        line = f"{label:<{max_name_len}}" 
        for res in residues:
            line += f'{res:<7}'
        score_file.write(line.rstrip() + '\n')


# --- Main Execution ---

def main():
    
    # --- Set Parameter Defaults ---
    group_id_cutoff = None 
    group_definition = None # Stores argument for -k (int or filename)
    target_k_groups = 0    
    group_gap_cutoff = 0.3
    column_gap_cutoff = 0.1
    cons_win_len = 3
    lamb = 0.7
    map_scores_to_01 = False
    matrix_file = ''
    outfile_prefix = ''
    align_file = None

    # --- Parse Options and Args ---
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hno:m:w:l:c:g:t:k:")
    except getopt.GetoptError as err:
        print(f"Error parsing arguments: {err}", file=sys.stderr); usage(); sys.exit(1)

    if len(args) != 1: usage(); sys.exit(2)
        
    align_file = args[0]
    outfile_prefix = align_file.split('.')[0]

    for opt, arg in opts:
        if opt == "-h": usage(); sys.exit()
        elif opt == "-n": map_scores_to_01 = True
        elif opt == "-o": outfile_prefix = arg
        elif opt == "-m": matrix_file = arg
        elif opt == "-w": 
            try: cons_win_len = int(arg)
            except ValueError: pass
        elif opt == "-l":
            try:
                if not (0.0 <= float(arg) <= 1.0): raise ValueError
                lamb = float(arg)
            except ValueError: pass
        elif opt == "-c":
            try:
                if not (0.0 <= float(arg) < 1.0): raise ValueError
                column_gap_cutoff = float(arg)
            except ValueError: pass
        elif opt == "-g":
            try:
                if not (0.0 <= float(arg) < 1.0): raise ValueError
                group_gap_cutoff = float(arg)
            except ValueError: pass
        elif opt == "-t":
            try:
                if not (0.0 <= float(arg) <= 100.0): raise ValueError
                group_id_cutoff = float(arg)
            except ValueError: pass
        elif opt == "-k":
            group_definition = arg
    
    # --- Validation and Group Mode Determination ---
    is_manual_k = group_definition and not group_definition.isdigit() # Check if -k is a filename
    is_auto_k = False
    
    if group_definition and group_definition.isdigit():
        try:
            target_k_groups = int(group_definition)
            if target_k_groups >= 2: is_auto_k = True
        except ValueError:
            pass # Should not happen if isdigit() is true

    if group_id_cutoff is None and not is_manual_k and not is_auto_k:
        print("\nERROR: You must specify EITHER an identity cutoff (-t), a target number of groups (-k N), OR a manual group file (-k filename).\n", file=sys.stderr)
        usage(); sys.exit(2)

    # --- Read Matrix and Alignment ---
    matrix, matrix_name = read_scoring_matrix(matrix_file)
    names, alignment = read_alignment(align_file)
    
    # --- Automated/Manual Group Inference ---
    similarity_matrix_df, linkage_matrix, final_identity_cutoff = None, None, 0.0

    if is_manual_k:
        print("\nStarting analysis using MANUALLY DEFINED groups from file...", file=sys.stderr)
        gids_to_seqs, final_identity_cutoff = parse_manual_groups_from_file(group_definition, names)
        
    else:
        print(f"\nStarting hierarchical clustering...", file=sys.stderr)
        similarity_matrix_df = calculate_similarity_matrix(names, alignment)
        dist_matrix = 1.0 - (similarity_matrix_df.values / 100.0)
        linkage_matrix = hierarchy.linkage(squareform(dist_matrix), method='average')
        
        if is_auto_k:
            # K-based clustering
            clusters = hierarchy.fcluster(linkage_matrix, target_k_groups, criterion='maxclust')
            cut_distance = linkage_matrix[len(names) - target_k_groups, 2]
            final_identity_cutoff = 100.0 * (1.0 - cut_distance)
            print(f"Target K={target_k_groups} groups achieved. Inferred identity cutoff: {final_identity_cutoff:.1f}%.", file=sys.stderr)
            
        else: 
            # T-based clustering
            group_id_cutoff = group_id_cutoff if group_id_cutoff is not None else 50.0
            cut_distance = 1.0 - (group_id_cutoff / 100.0)
            clusters = hierarchy.fcluster(linkage_matrix, cut_distance, criterion='distance')
            final_identity_cutoff = group_id_cutoff

        # Map clusters to groups
        gids_to_seqs = {}
        for i, cluster_id in enumerate(clusters):
            gid = f"Group_{cluster_id}"
            if gid not in gids_to_seqs: gids_to_seqs[gid] = []
            gids_to_seqs[gid].append(i)

    group_ids = sorted(gids_to_seqs.keys())
    num_groups = len(group_ids)

    if num_groups < 2:
        sys.exit(f"ERROR: Grouping resulted in only {num_groups} group(s). Please check your input parameters.")

    # --- Score Columns with GroupSim ---
    scores = []
    num_cols = len(alignment[0])

    for i in range(num_cols):
        col = get_column(i, alignment)
        groups_of_residues = [[col[seq_num] for seq_num in gids_to_seqs[gid]] for gid in group_ids]
        group_too_gappy = any(gap_percentage(g) > group_gap_cutoff for g in groups_of_residues)

        if group_too_gappy or gap_percentage(col) > column_gap_cutoff:
            col_score = None
        else:
            col_score = matrix_sum_pairs(groups_of_residues, matrix)
        scores.append(col_score)

    if map_scores_to_01: scores = norm_01(scores)
    if cons_win_len > 0: scores = cons_window_score(scores, alignment, cons_win_len, lamb)


    # --- Visualization ---
    if linkage_matrix is not None:
        plot_clustered_heatmap(similarity_matrix_df, linkage_matrix, outfile_prefix)
        plot_dendrogram(linkage_matrix, names, final_identity_cutoff, outfile_prefix)

    plot_groupsim_manhattan(scores, alignment, gids_to_seqs, outfile_prefix)


    # --- Print Scores to Output ---
    score_filename = f"{outfile_prefix}.txt"
    try:
        score_file = open(score_filename, 'w')
    except IOError:
        sys.exit(f"Error: Could not open output file {score_filename}.")

    # Write headers
    grouping_method = "Manual" if is_manual_k else "Clustering"
    score_file.write(f'# {align_file} scored by GroupSim with {matrix_name} matrix (Grouping by: {grouping_method})\n')
    score_file.write(f'# Grouping criterion: {final_identity_cutoff:.1f}% identity cutoff\n')
    score_file.write(f'# window params: len={cons_win_len}, lambda={lamb:.2f}\n')
    score_file.write(f'# group gap cutoff: {group_gap_cutoff:.2f}  column gap cutoff: {column_gap_cutoff:.2f}\n')
    score_file.write(f'# Groups: {", ".join(group_ids)}\n')
    score_file.write('# col_num\tscore\tcolumn_residues\n')

    # Write column scores
    for i, col_score in enumerate(scores):
        col_str = format_column(i, alignment, gids_to_seqs)
        score_file.write(f'{i}\t{col_score:.6f}\t{col_str}\n' if col_score is not None else f'{i}\tNone\t{col_str}\n')

    # Detailed Output
    write_detailed_output(score_file, names, alignment, scores, gids_to_seqs, group_ids, final_identity_cutoff)

    score_file.close()
    print(f"\nWrote scores and detailed output to {score_filename}.", file=sys.stderr)


if __name__ == "__main__":
    main()
# End of auto_groupsim_final_v3.py
