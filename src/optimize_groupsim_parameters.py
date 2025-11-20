#!/usr/bin/env python
"""
optimize_groupsim_params.py - Demonstrates a grid search for optimal GroupSim parameters.
Analyzes the Z-scores of known SDPs (Pos 15) to find the best (w, l) combination.
"""
from scipy.stats import zscore 
import subprocess
import numpy as np
import pandas as pd
import re
import sys

# Define the ranges to test
W_RANGES = [1, 3, 5, 7]
L_RANGES = [0.5, 0.6, 0.7, 0.8, 0.9]
TARGET_SDP_POSITIONS = [15]  # Known perfect SDP from synthetic data

MAIN_SCRIPT = "groupsim.py"
ALIGNMENT_FILE = "synthetic_test.fasta"
GROUP_FILE = "synthetic_test_groups.txt"
OUTPUT_PREFIX = "optimization_run"

def run_groupsim(w, l, align_file, group_file, output_prefix):
    """Executes the main GroupSim script and returns the score file content."""
    
    cmd = [
        "python", MAIN_SCRIPT, 
        "-k", group_file, 
        "-w", str(w), 
        "-l", str(l), 
        "-o", f"{output_prefix}_w{w}_l{l:.1f}",
        align_file
    ]
    
    try:
        # Run the command
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        score_filename = f"{output_prefix}_w{w}_l{l:.1f}.txt"
        
        # Read the scores from the output file
        with open(score_filename, 'r') as f:
            content = f.read()
        return content

    except subprocess.CalledProcessError as e:
        print(f"Error running GroupSim for w={w}, l={l}: {e.stderr}", file=sys.stderr)
        return None
    except FileNotFoundError:
        print(f"Error: Could not find {MAIN_SCRIPT}. Ensure it is in the current directory.", file=sys.stderr)
        sys.exit(1)


def parse_z_scores(score_file_content, target_positions):
    """Parses scores, calculates Z-scores, and extracts the Z-score for target positions."""
    
    scores_list = []
    # Find the line containing scores (starts with '# col_num')
    score_data_pattern = re.compile(r'^\d+\t(.*?)\t.*$', re.MULTILINE)
    
    for line in score_file_content.split('\n'):
        if line.startswith('#') or not line.strip():
            continue
        try:
            parts = line.split('\t')
            pos = int(parts[0])
            score_str = parts[1].strip()
            
            if score_str != 'None':
                scores_list.append((pos, float(score_str)))
        except IndexError:
            continue
        except ValueError:
            continue

    if not scores_list:
        return {pos: 0.0 for pos in target_positions}, 0.0

    scores_df = pd.DataFrame(scores_list, columns=['pos', 'score'])
    
    # Calculate global Z-scores
    scores_df['z_score'] = zscore(scores_df['score'])
    
    # Extract Z-score for target positions
    target_z_scores = {}
    for pos in target_positions:
        z_score = scores_df[scores_df['pos'] == pos]['z_score'].iloc[0] if not scores_df[scores_df['pos'] == pos].empty else 0.0
        target_z_scores[pos] = z_score
        
    # Metric for optimization: Average Z-score of target SDPs
    avg_target_z = np.mean(list(target_z_scores.values())) if target_z_scores else 0.0
    
    return target_z_scores, avg_target_z


def main_optimization():
    
    if not all([os.path.exists(MAIN_SCRIPT), os.path.exists(ALIGNMENT_FILE), os.path.exists(GROUP_FILE)]):
        print("Please ensure you have run 'generate_synthetic_alignment.py' and that 'groupsim.py' is in your directory.")
        return

    results = []
    
    for w in W_RANGES:
        for l in L_RANGES:
            print(f"--- Running w={w}, l={l:.1f} ---", file=sys.stderr)
            
            content = run_groupsim(w, l, ALIGNMENT_FILE, GROUP_FILE, OUTPUT_PREFIX)
            
            if content:
                target_z_scores, avg_target_z = parse_z_scores(content, TARGET_SDP_POSITIONS)
                
                # Add individual SDP Z-scores to the result list
                result = {'W': w, 'L': l, 'Avg_Z_SDPs': avg_target_z}
                for pos, z in target_z_scores.items():
                    result[f'Z_Pos_{pos}'] = z
                results.append(result)

    if results:
        df = pd.DataFrame(results)
        df_sorted = df.sort_values(by='Avg_Z_SDPs', ascending=False)
        
        print("\n" + "="*80)
        print("                 GroupSim Parameter Optimization Results (Grid Search)")
        print("="*80)
        print(df_sorted.to_markdown(index=False, floatfmt=".3f"))
        print("\nRecommendation: Choose the (W, L) pair that maximizes Avg_Z_SDPs.")
    else:
        print("Optimization failed to produce results.")

import os
if __name__ == "__main__":
    if not os.path.exists(ALIGNMENT_FILE) or not os.path.exists(GROUP_FILE):
        print(f"Missing {ALIGNMENT_FILE} or {GROUP_FILE}. Please run generate_synthetic_alignment.py first.")
    main_optimization()
