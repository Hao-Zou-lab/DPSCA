import argparse
import os
import pandas as pd
from Bio import Phylo
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Find best gene matches based on Phylogeny, Orthology, and BLAST.")
    parser.add_argument("--tree_dir", required=True, help="Directory containing tree files (e.g., OG0000010_tree.txt)")
    parser.add_argument("--og_file", required=True, help="Orthogroup table (Species columns, OG rows)")
    parser.add_argument("--blast_file", required=True, help="BLAST result file (outfmt 6)")
    parser.add_argument("--target_ids", required=True, help="File with list of Target Gene IDs")
    parser.add_argument("--target_species", required=True, help="Name of the target species (as used in OG header and Tree prefix)")
    parser.add_argument("--ref_species", required=True, help="Comma-separated list of reference species names")
    parser.add_argument("--output", required=True, help="Output file path")
    return parser.parse_args()

def load_og_data(og_file, target_species, ref_species_list):
    """
    Parses the OG table.
    Returns:
    1. gene_to_og: {gene_id: og_id}
    2. og_to_species_genes: {og_id: {species: [gene_list]}}
    3. gene_to_species: {gene_id: species} (helper to identify species of blast hits)
    """
    print("Loading OG table...")
    df = pd.read_csv(og_file, sep='\t', dtype=str)
    
    gene_to_og = {}
    og_to_species_genes = {}
    gene_to_species = {}
    
    # Check headers
    headers = df.columns.tolist()
    # Assuming first column is OG ID
    og_col = headers[0]
    
    valid_species = [target_species] + ref_species_list
    
    for _, row in df.iterrows():
        og_id = row[og_col]
        if pd.isna(og_id): continue
        
        og_to_species_genes[og_id] = {}
        
        for sp in valid_species:
            if sp not in headers:
                continue
                
            genes_str = row[sp]
            if pd.isna(genes_str):
                og_to_species_genes[og_id][sp] = []
                continue
            
            # Split by comma and strip whitespace
            genes = [g.strip() for g in genes_str.split(',')]
            og_to_species_genes[og_id][sp] = genes
            
            for g in genes:
                gene_to_og[g] = og_id
                gene_to_species[g] = sp
                
    return gene_to_og, og_to_species_genes, gene_to_species

def load_blast(blast_file):
    """
    Loads BLAST output (fmt 6). 
    Returns: {query_gene: {subject_gene: bitscore}}
    """
    print("Loading BLAST results...")
    blast_scores = {}
    try:
        # Read only necessary columns: qseqid (0), sseqid (1), bitscore (11)
        # Using pandas for speed on large files
        df = pd.read_csv(blast_file, sep='\t', header=None, usecols=[0, 1, 11], names=['q', 's', 'score'])
        
        # Group by query to create dictionary
        for q_gene, group in df.groupby('q'):
            blast_scores[q_gene] = dict(zip(group['s'], group['score']))
            
    except Exception as e:
        print(f"Error reading BLAST file: {e}")
        sys.exit(1)
        
    return blast_scores

def get_tree_distances(tree_path, target_species, target_gene, candidate_genes, ref_species):
    """
    Parses the tree and calculates distance from target to all candidates.
    Returns a dict: {ref_gene: distance}
    Note: Tree nodes are "Species_Gene", inputs are "Gene".
    """
    if not os.path.exists(tree_path):
        return {}

    try:
        tree = Phylo.read(tree_path, "newick")
    except:
        return {}

    # Helper to find node name in tree (handling the Species_ prefix)
    def find_node(species, gene):
        target_name = f"{species}_{gene}"
        # Exact match search
        matches = [n for n in tree.get_terminals() if n.name == target_name]
        if matches:
            return matches[0]
        # Fallback: sometimes tree names might vary slightly, but assuming strict format per prompt
        return None

    target_node = find_node(target_species, target_gene)
    if not target_node:
        return {}

    distances = {}
    for gene in candidate_genes:
        # candidate gene belongs to ref_species
        ref_node = find_node(ref_species, gene)
        if ref_node:
            try:
                # distance method calculates branch length distance
                dist = tree.distance(target_node, ref_node)
                distances[gene] = dist
            except:
                pass
    
    return distances

def main():
    args = parse_args()
    
    ref_species_list = args.ref_species.split(',')
    
    # 1. Load Data
    gene_to_og, og_to_species_genes, gene_to_species = load_og_data(args.og_file, args.target_species, ref_species_list)
    blast_scores = load_blast(args.blast_file)
    
    with open(args.target_ids, 'r') as f:
        target_genes = [line.strip() for line in f if line.strip()]

    # 2. Process
    print("Processing target genes...")
    results = []
    
    # Header for output
    header = ["Target_Gene", "OG_ID"]
    for ref in ref_species_list:
        header.extend([f"{ref}_Best_Match", f"{ref}_Score", f"{ref}_Other_Candidates"])
    
    results.append("\t".join(header))

    for t_gene in target_genes:
        row = [t_gene]
        
        og_id = gene_to_og.get(t_gene, "NA")
        row.append(og_id)
        
        # Determine tree file path if OG exists
        tree_path = None
        if og_id != "NA":
            tree_path = os.path.join(args.tree_dir, f"{og_id}_tree.txt")
        
        # Get BLAST hits for this target
        t_blast_hits = blast_scores.get(t_gene, {}) # {subject_gene: score}
        
        for ref_sp in ref_species_list:
            # Step A: Identify Candidates
            # 1. From OG
            og_candidates = set()
            if og_id != "NA" and og_id in og_to_species_genes:
                og_candidates = set(og_to_species_genes[og_id].get(ref_sp, []))
            
            # 2. From BLAST (only if they belong to this ref_species)
            blast_candidates = set()
            for hit_gene in t_blast_hits.keys():
                # We need to know if this hit_gene belongs to ref_sp
                # We use the gene_to_species map derived from OG file
                if gene_to_species.get(hit_gene) == ref_sp:
                    blast_candidates.add(hit_gene)
            
            # Union of candidates
            all_candidates = list(og_candidates.union(blast_candidates))
            
            if not all_candidates:
                row.extend(["NA", "0", "NA"])
                continue
            
            # Step B: Pre-calculate metrics for scoring
            
            # 1. Tree Distances
            # Only if we have a tree and candidates
            tree_dists = {}
            if tree_path:
                tree_dists = get_tree_distances(tree_path, args.target_species, t_gene, all_candidates, ref_sp)
            
            min_dist = float('inf')
            if tree_dists:
                min_dist = min(tree_dists.values())

            # 2. Max BLAST Score
            # specific to this ref species
            current_blast_scores = {g: t_blast_hits.get(g, 0) for g in all_candidates}
            max_bitscore = 0
            if current_blast_scores:
                max_bitscore = max(current_blast_scores.values())

            # 3. OG Count (for single copy check)
            # Only counts genes actually in the OG, not extra BLAST hits
            og_count = len(og_candidates)

            # Step C: Calculate Scores
            scored_candidates = []
            
            for cand in all_candidates:
                score = 0.0
                
                # Rule 1: Same OG (0.2)
                if cand in og_candidates:
                    score += 0.2
                    
                # Rule 2: Phylogenetically Closest (0.4)
                # Must be in tree and have distance equal to min_dist
                if cand in tree_dists and tree_dists[cand] == min_dist:
                    score += 0.4
                    
                # Rule 3: Highest BLAST Score (0.35)
                # Must have score > 0 and equal to max
                b_score = current_blast_scores.get(cand, 0)
                if b_score > 0 and b_score == max_bitscore:
                    score += 0.35
                    
                # Rule 4: Single Copy in Ref (0.05)
                # Logic: If the OG has exactly 1 gene for this species, that gene gets points
                if og_count == 1 and cand in og_candidates:
                    score += 0.05
                
                scored_candidates.append({
                    "gene": cand,
                    "score": round(score, 3),
                    "blast": b_score,
                    "dist": tree_dists.get(cand, float('inf'))
                })
            
            # Step D: Select Best Match
            # Sort by Score (desc), then BLAST (desc), then Distance (asc)
            scored_candidates.sort(key=lambda x: (-x['score'], -x['blast'], x['dist']))
            
            if not scored_candidates:
                row.extend(["NA", "0", "NA"])
            else:
                best = scored_candidates[0]
                others = scored_candidates[1:]
                
                # Format Best: "GeneID(Score:X)"
                best_str = f"{best['gene']}"
                best_score_str = str(best['score'])
                
                # Format Others: "GeneID(Score:X);..."
                others_str = ";".join([f"{o['gene']}({o['score']})" for o in others])
                if not others_str:
                    others_str = "None"
                
                row.extend([best_str, best_score_str, others_str])

        results.append("\t".join(row))

    # Write Output
    with open(args.output, 'w') as f:
        f.write("\n".join(results))
    
    print(f"Done! Results written to {args.output}")

if __name__ == "__main__":
    main()
