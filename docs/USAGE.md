# GroupSim-Py3 Usage Guide

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Usage Modes](#usage-modes)
4. [Parameter Guide](#parameter-guide)
5. [Output Files](#output-files)
6. [Workflow Examples](#workflow-examples)

## Installation

### 1. Clone Repository

```bash
git clone git@github.com:jacgonisa/groupsim-py3.git
cd groupsim-py3
```

### 2. Create Environment

```bash
conda create -n groupsim python=3.8
conda activate groupsim
pip install -r requirements.txt
```

### 3. Test Installation

```bash
cd data/examples
python ../../src/groupsim.py -k synthetic_test_groups.txt synthetic_test.fasta
```

Expected output:
- `synthetic_test_manhattan_plot.png`
- `synthetic_test.txt`

## Quick Start

### Scenario 1: I have an alignment, want automatic grouping

```bash
# Group sequences with ≥75% identity
python groupsim.py -t 75.0 my_alignment.fasta
```

### Scenario 2: I know I want exactly 3 groups

```bash
# Algorithm finds the cutoff for 3 groups
python groupsim.py -k 3 my_alignment.fasta
```

### Scenario 3: I have predefined groups

Create `my_groups.txt`:
```
Family_A: Seq1, Seq2, Seq3
Family_B: Seq4, Seq5, Seq6
```

Run:
```bash
python groupsim.py -k my_groups.txt my_alignment.fasta
```

## Usage Modes

### Mode 1: Identity Cutoff (`-t`)

**Use when**: You know the similarity threshold that defines your groups.

```bash
python groupsim.py -t 80.0 alignment.fasta
```

- Groups sequences with ≥80% pairwise identity
- Produces dendrogram showing the cutoff line
- Best for well-separated functional groups

**Choosing a cutoff**:
- 90-100%: Very closely related (species/strains)
- 70-90%: Ortholog groups
- 50-70%: Paralog families
- <50%: Distantly related subfamilies

### Mode 2: Target K Groups (`-k N`)

**Use when**: You know the number of functional groups but not the threshold.

```bash
python groupsim.py -k 4 alignment.fasta
```

- Creates exactly 4 groups
- Algorithm determines optimal identity cutoff
- Good for exploratory analysis

### Mode 3: Manual Groups (`-k file`)

**Use when**: Groups are defined by prior knowledge (e.g., experimental data, literature).

**Group file format**:
```
# Comments start with #
Group_Kinase: Seq1, Seq2, Seq5
Group_Phosphatase: Seq3, Seq4
# Unlisted sequences → 'Other' group automatically
```

```bash
python groupsim.py -k groups.txt alignment.fasta
```

## Parameter Guide

### Essential Parameters

| Parameter | Description | Typical Values |
|-----------|-------------|----------------|
| `-t` | Identity cutoff (%) | 50-90 |
| `-k` | Target groups or file | 2-10 groups |
| `-o` | Output prefix | any string |

### Scoring Parameters

| Parameter | Description | Default | When to Change |
|-----------|-------------|---------|----------------|
| `-w` | Window size | 3 | Increase (5-7) for long-range effects |
| `-l` | Lambda (context weight) | 0.7 | Decrease (0.5) if context less important |
| `-c` | Column gap cutoff | 0.1 | Increase (0.3) for gappy alignments |
| `-g` | Group gap cutoff | 0.3 | Increase (0.5) for highly variable groups |

### Advanced Options

| Parameter | Description | Use Case |
|-----------|-------------|----------|
| `-m` | Similarity matrix file | Use BLOSUM62 for biochemical similarity |
| `-n` | Normalize scores [0,1] | For comparing across alignments |

## Output Files

### 1. Manhattan Plot (`*_manhattan_plot.png`)

**What it shows**:
- X-axis: Alignment position
- Y-axis: GroupSim score
- Color: Z-score (red = high specificity)
- Annotations: Group-specific residues at significant positions

**How to read**:
- **High peaks**: Strong SDPs
- **Red points**: Z-score > threshold (default 2.0)
- **Annotations**: Show residue pattern (e.g., "G1:K | G2:D")

### 2. Clustered Heatmap (`*_heatmap.png`)

**What it shows**:
- Pairwise sequence similarities
- Hierarchical clustering (dendrogram on sides)
- Group structure

**How to read**:
- **Yellow**: High similarity
- **Purple**: Low similarity
- **Blocks**: Natural groupings
- **Dendrogram height**: Divergence between clusters

### 3. Dendrogram (`*_dendrogram.png`)

**What it shows**:
- Hierarchical tree of sequences
- Red dashed line: Identity cutoff used for grouping
- Top axis: Percent identity
- Bottom axis: Distance (1 - identity/100)

**How to read**:
- **Clusters below line**: Belong to same group
- **Branch length**: Evolutionary/sequence distance
- **Leaf colors**: Indicate group membership

### 4. Score File (`*.txt`)

**Contents**:
```
# Header: Parameters used
# Column scores: Position, score, residues per group
# Detailed alignment: Full alignment with group labels
```

**Key sections**:
1. **Summary**: Parameters and group definitions
2. **Column scores**: One line per position
3. **Detailed alignment**: Visual alignment with scores

## Workflow Examples

### Example 1: Kinase Subfamily Analysis

**Goal**: Find residues determining substrate specificity in kinase subfamilies.

```bash
# Step 1: Align sequences (outside GroupSim)
mafft --auto kinases.fasta > kinases_aligned.fasta

# Step 2: Run GroupSim (assuming 3 subfamilies)
python groupsim.py -k 3 kinases_aligned.fasta -o kinase_analysis

# Step 3: Extract top SDPs
awk '$2 != "None" && $2 > 0.5' kinase_analysis.txt | sort -k2 -rn | head -20

# Step 4: Review Manhattan plot
# → Look for positions near active site or substrate-binding pocket

# Step 5: Validate experimentally
# → Mutagenesis of top SDP positions
```

### Example 2: Optimization for Your Data

**Goal**: Find best parameters for your specific alignment.

```bash
# Step 1: Identify known SDPs (from literature, structure, etc.)
# Example: Positions 42, 108, 215

# Step 2: Edit optimizer
nano src/optimize_groupsim_parameters.py
# Change:
# TARGET_SDP_POSITIONS = [42, 108, 215]
# ALIGNMENT_FILE = "your_alignment.fasta"
# GROUP_FILE = "your_groups.txt"

# Step 3: Run optimization
cd src/
python optimize_groupsim_parameters.py

# Step 4: Check results
# → Table shows Z-scores for each (W, L) combination
# → Choose pair with highest Avg_Z_SDPs

# Step 5: Run with optimal parameters
python groupsim.py -k your_groups.txt -w 5 -l 0.8 your_alignment.fasta
```

### Example 3: Large-Scale Family Analysis

**Goal**: Analyze 50+ sequences, multiple subfamilies.

```bash
# Step 1: Prepare groups (if manual mode)
# Create groups.txt with subfamily definitions

# Step 2: Run with adjusted gap thresholds
python groupsim.py -k groups.txt -c 0.2 -g 0.4 large_family.fasta -o family_sdp

# Higher gap thresholds handle diverse alignments

# Step 3: Filter high-confidence SDPs
awk '$2 != "None"' family_sdp.txt | \
  awk '{z = ($2 - mean) / sd} z > 3.0' > top_sdps.txt
# (Calculate mean and sd from scores first)

# Step 4: Visualize in protein structure
# Map SDP positions to 3D structure using PyMOL:
# select sdp_sites, resi 42+108+215
# show spheres, sdp_sites
# color red, sdp_sites
```

## Common Issues and Solutions

### Issue: All scores are very low

**Likely cause**: Groups are too similar or parameters too strict.

**Solutions**:
1. Check if groups are actually different (view heatmap)
2. Decrease `-c` and `-g` (allow more gaps)
3. Increase `-w` (larger context window)
4. Try `-n` flag (normalize scores)

### Issue: Only a few positions score high

**This may be correct!** True SDPs can be rare. However, if unexpected:

**Solutions**:
1. Adjust `-l` (lambda): Try 0.5-0.9 range
2. Check alignment quality (realign if needed)
3. Verify group definitions are biologically meaningful

### Issue: Too many high-scoring positions

**Likely cause**: Groups are very divergent or lambda too low.

**Solutions**:
1. Increase `-l` to weight raw scores more (e.g., 0.9)
2. Check if groups should be subdivided
3. Use stricter gap thresholds (`-c 0.05 -g 0.2`)

## Tips for Best Results

### 1. Alignment Quality Matters
- Remove highly divergent/fragmentary sequences before alignment
- Use profile alignment methods (e.g., MAFFT --auto) for diverse families
- Manually inspect alignment in critical regions

### 2. Group Definition Strategy
- **Start automatic**: Use `-k N` to explore natural groupings
- **Refine manually**: Based on phylogeny, function, or structure
- **Validate biologically**: Do groups correspond to known subfamilies?

### 3. Parameter Optimization
- Default parameters (w=3, l=0.7) work well for most cases
- Optimize only if you have known SDPs to validate against
- Document parameter choices for reproducibility

### 4. Interpreting Scores
- **Don't rely on scores alone**: Consider structural/functional context
- **Validate experimentally**: Top SDPs are predictions, not proof
- **Look for clusters**: Neighboring high-scoring positions often form functional sites

## Next Steps

After identifying SDPs:

1. **Structural mapping**: Visualize SDPs on 3D structure
2. **Functional validation**: Site-directed mutagenesis
3. **Co-evolution analysis**: Check for correlated mutations
4. **Predictive modeling**: Use SDPs for subfamily classification
5. **Literature mining**: Check if SDPs are known specificity determinants

## Support

For issues, questions, or suggestions:
- **GitHub Issues**: https://github.com/jacgonisa/groupsim-py3/issues
- **Documentation**: See README.md for detailed methodology
- **Original paper**: Capra & Singh (2008) for theoretical background
