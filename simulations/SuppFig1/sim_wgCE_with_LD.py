#!/usr/bin/env python3
"""
Simulation: wg-CE calibration with realistic LD — clumped-PRS design.

Design:
  - Load per-chromosome non-LD-pruned .bed files (real LD structure preserved)
  - Select M causal SNPs randomly across all chromosomes
    (causal SNPs may be in LD with each other — this is realistic)
  - Select P total SNPs (default 100,000) including all M causal SNPs
    (simulates an LD-clumped PRS: mostly independent, but not perfectly)
  - Load N_gwas + N_test non-overlapping individuals
  - Simulate phenotype using only the M causal SNPs (two-pathway model)
  - Run GWAS on ALL P SNPs in the GWAS sample (N_gwas individuals)
  - Build PRS on ALL P SNPs in the independent test sample (N_test individuals)
  - Run wg-CE and pc-CE on the test sample
  - Compare to oracle CE (using true liabilities)

Usage:
    python sim_wgCE_with_LD.py \
        --bed-prefix '/path/to/ukb_imp_v3_chr{CHR}.qced' \
        --chromosomes 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
        --out results_LD_clumped --n-causal 1000 --n-prs-snps 100000 \
        --n-gwas 10000 --n-test 10000 --n-sim 500

    Expects files like: {prefix}1.bed/.bim/.fam or with {CHR} placeholder

Requirements:
    pip install numpy pandas scipy matplotlib bed-reader statsmodels
"""

import argparse
import os
import numpy as np
import pandas as pd
from scipy import stats
import warnings
import time

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# 1. Genotype loading — random genome-wide sampling (no LD-window logic)
# ---------------------------------------------------------------------------

def _find_chr_bed(prefix, chrom):
    """
    Find the .bed file for a given chromosome.

    If prefix contains '{CHR}', replace it with the chromosome number.
      e.g. '/path/to/ukb_imp_v3_chr{CHR}.qced.bed' -> '.../ukb_imp_v3_chr10.qced.bed'

    Otherwise try common suffix patterns:
      {prefix}{chrom}.bed, {prefix}chr{chrom}.bed, {prefix}.chr{chrom}.bed
    """
    if "{CHR}" in prefix:
        path = prefix.replace("{CHR}", str(chrom))
        if not path.endswith(".bed"):
            path += ".bed"
        if os.path.exists(path):
            return path
        return None

    candidates = [
        f"{prefix}{chrom}.bed",
        f"{prefix}chr{chrom}.bed",
        f"{prefix}.chr{chrom}.bed",
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    return None


def load_genotypes(bed_path=None, bed_prefix=None, chromosomes=None,
                   n_samples=20000, n_causal=1000, n_prs_snps=100000,
                   maf_min=0.05, seed=42):
    """
    Load real genotypes for clumped-PRS simulation.

    Strategy (per-chromosome mode):
      1. Read .bim files, screen for MAF to get common SNPs
      2. Pick n_causal random common SNPs as causal (across all chromosomes)
      3. Pick n_prs_snps total SNPs (including all causal) for PRS construction
      4. Load genotypes for the n_prs_snps SNPs

    Returns:
        G_all      : (n_samples, n_prs_snps) standardised genotype matrix
        G_causal   : (n_samples, n_causal)   causal genotype matrix
        causal_mask: boolean array of length n_prs_snps
        snp_info   : DataFrame with chr, pos, is_causal for ALL PRS SNPs
        causal_info: DataFrame for causal SNPs only
        sample_idx : indices of selected individuals
    """
    rng = np.random.default_rng(seed)

    if bed_path is not None:
        return _load_single(bed_path, n_samples, n_causal, n_prs_snps,
                            maf_min, rng)
    elif bed_prefix is not None:
        if chromosomes is None:
            chromosomes = [str(c) for c in range(1, 23)]
        return _load_perchr(bed_prefix, chromosomes, n_samples,
                            n_causal, n_prs_snps, maf_min, rng)
    else:
        raise ValueError("Provide either --bed or --bed-prefix")


def _load_single(bed_path, n_samples, n_causal, n_prs_snps, maf_min, rng):
    """Load from a single .bed file."""
    from bed_reader import open_bed

    # Read .bim for positions
    bim_path = str(bed_path).replace(".bed", ".bim")
    if not os.path.exists(bim_path):
        raise FileNotFoundError(f"Need .bim file: {bim_path}")
    bim = pd.read_csv(bim_path, sep="\t", header=None,
                      names=["chr", "snp", "cm", "pos", "a1", "a2"],
                      dtype={"chr": str})

    print(f"  Opening {bed_path}...")
    with open_bed(bed_path) as bed:
        n_tot_samples = bed.iid_count
        n_tot_snps = bed.sid_count
        print(f"  BED: {n_tot_samples} samples, {n_tot_snps} SNPs")

        if n_tot_samples < n_samples:
            raise ValueError(f"Need {n_samples} samples, only {n_tot_samples}")

        sample_idx = np.sort(
            rng.choice(n_tot_samples, size=n_samples, replace=False)
        ).astype(np.intp)

        # Screen ALL SNPs for MAF (in batches)
        print(f"  Screening all {n_tot_snps} SNPs for MAF...")
        batch_size = 50000
        maf_arr = np.zeros(n_tot_snps)
        for start in range(0, n_tot_snps, batch_size):
            end = min(start + batch_size, n_tot_snps)
            idx_batch = list(range(start, end))
            g = bed.read(index=(sample_idx.tolist(), idx_batch), dtype='float64')
            g = np.nan_to_num(g, nan=0.0)
            af = g.mean(axis=0) / 2.0
            maf_arr[start:end] = np.minimum(af, 1.0 - af)

        common_indices = np.where(maf_arr >= maf_min)[0]
        n_common = len(common_indices)
        print(f"  {n_common} common SNPs (MAF >= {maf_min})")

        if n_common < n_prs_snps:
            raise ValueError(f"Need {n_prs_snps} PRS SNPs but only {n_common} common SNPs")

        # Pick causal SNPs from common ones
        causal_in_common = np.sort(rng.choice(n_common, size=n_causal, replace=False))
        causal_global = common_indices[causal_in_common]

        # Pick remaining PRS SNPs (excluding causal, from common)
        n_extra = n_prs_snps - n_causal
        non_causal_common = np.setdiff1d(np.arange(n_common), causal_in_common)
        extra_in_common = np.sort(rng.choice(len(non_causal_common), size=n_extra, replace=False))
        extra_global = common_indices[non_causal_common[extra_in_common]]

        # Combined PRS SNP set
        all_global = np.sort(np.concatenate([causal_global, extra_global]))
        causal_mask = np.isin(all_global, causal_global)

        print(f"  Selected {len(all_global)} PRS SNPs ({n_causal} causal + {n_extra} non-causal)")

        # Read genotypes
        print(f"  Reading genotypes for {len(all_global)} SNPs x {n_samples} individuals...")
        G_all = bed.read(
            index=(sample_idx.tolist(), all_global.tolist()), dtype='float64'
        )
        G_all = np.nan_to_num(G_all, nan=0.0)

    snp_chrs = bim.iloc[all_global]["chr"].values
    snp_pos = bim.iloc[all_global]["pos"].values

    return _finalise(G_all, snp_chrs, snp_pos, causal_mask, sample_idx)


def _load_perchr(bed_prefix, chromosomes, n_samples, n_causal, n_prs_snps,
                 maf_min, rng):
    """
    Load from per-chromosome .bed files.

      1. Read all .bim files, screen for MAF
      2. Pick n_causal random common SNPs as causal
      3. Pick n_prs_snps total (including causal) for PRS
      4. Load genotypes chromosome by chromosome
    """
    from bed_reader import open_bed

    # --- Phase 1: discover files, read .bim files ---
    chr_files = {}
    chr_bims = {}
    for chrom in chromosomes:
        path = _find_chr_bed(bed_prefix, chrom)
        if path is None:
            print(f"  WARNING: no bed file for chr {chrom}, skipping")
            continue
        chr_files[chrom] = path
        bim_path = path.replace(".bed", ".bim")
        if not os.path.exists(bim_path):
            raise FileNotFoundError(f"Need .bim file: {bim_path}")
        chr_bims[chrom] = pd.read_csv(
            bim_path, sep="\t", header=None,
            names=["chr", "snp", "cm", "pos", "a1", "a2"],
            dtype={"chr": str}
        )

    if not chr_files:
        raise FileNotFoundError(f"No .bed files found with prefix '{bed_prefix}'")

    sorted_chroms = sorted(chr_files.keys(), key=lambda c: int(c) if c.isdigit() else 99)
    print(f"  Found {len(chr_files)} chromosome files")

    # --- Phase 2: pick sample indices from first chromosome ---
    first_path = chr_files[sorted_chroms[0]]
    with open_bed(first_path) as bed:
        n_tot_samples = bed.iid_count
        if n_tot_samples < n_samples:
            raise ValueError(f"Need {n_samples} samples, only {n_tot_samples}")
        sample_idx = np.sort(
            rng.choice(n_tot_samples, size=n_samples, replace=False)
        ).astype(np.intp)
    print(f"  Selected {n_samples} individuals")

    # --- Phase 3: screen each chromosome for MAF ---
    chr_common_local = {}   # chrom -> array of local indices that are common
    chr_positions = {}      # chrom -> positions array (all SNPs)

    for chrom in sorted_chroms:
        bim_chr = chr_bims[chrom]
        n_snps_chr = len(bim_chr)

        with open_bed(chr_files[chrom]) as bed:
            batch_size = 50000
            maf_chr = np.zeros(n_snps_chr)
            for start in range(0, n_snps_chr, batch_size):
                end = min(start + batch_size, n_snps_chr)
                idx_batch = list(range(start, end))
                g = bed.read(index=(sample_idx.tolist(), idx_batch), dtype='float64')
                g = np.nan_to_num(g, nan=0.0)
                af = g.mean(axis=0) / 2.0
                maf_chr[start:end] = np.minimum(af, 1.0 - af)

        common_local = np.where(maf_chr >= maf_min)[0]
        chr_common_local[chrom] = common_local
        chr_positions[chrom] = bim_chr["pos"].values
        print(f"    Chr {chrom:>2}: {n_snps_chr:>7} total, {len(common_local):>6} common")

    total_common = sum(len(v) for v in chr_common_local.values())
    print(f"  Total common SNPs: {total_common}")

    if total_common < n_prs_snps:
        raise ValueError(f"Need {n_prs_snps} PRS SNPs but only {total_common} common SNPs")

    # --- Phase 4: pick causal SNPs (genome-wide, proportional to chr size) ---
    # Build a flat pool of (chrom, local_idx) for all common SNPs
    pool_chroms = []
    pool_local_idx = []
    for chrom in sorted_chroms:
        for idx in chr_common_local[chrom]:
            pool_chroms.append(chrom)
            pool_local_idx.append(idx)

    pool_chroms = np.array(pool_chroms)
    pool_local_idx = np.array(pool_local_idx)

    # Pick causal
    causal_pick = rng.choice(total_common, size=n_causal, replace=False)
    causal_chroms = pool_chroms[causal_pick]
    causal_local_indices = pool_local_idx[causal_pick]

    # Group causal by chromosome
    causal_by_chr = {}
    for i, chrom in enumerate(causal_chroms):
        causal_by_chr.setdefault(chrom, set()).add(causal_local_indices[i])

    n_causal_per_chr = {c: len(v) for c, v in sorted(causal_by_chr.items(),
                        key=lambda x: int(x[0]) if x[0].isdigit() else 99)}
    print(f"  Causal SNPs per chromosome: {n_causal_per_chr}")
    print(f"  Total causal: {sum(n_causal_per_chr.values())}")

    # --- Phase 5: pick remaining PRS SNPs ---
    # Remove causal from pool, then sample n_extra more
    n_extra = n_prs_snps - n_causal
    non_causal_pool = np.setdiff1d(np.arange(total_common), causal_pick)
    extra_pick = rng.choice(len(non_causal_pool), size=n_extra, replace=False)
    extra_pool_idx = non_causal_pool[extra_pick]
    extra_chroms = pool_chroms[extra_pool_idx]
    extra_local_indices = pool_local_idx[extra_pool_idx]

    # Group extra by chromosome
    extra_by_chr = {}
    for i, chrom in enumerate(extra_chroms):
        extra_by_chr.setdefault(chrom, set()).add(extra_local_indices[i])

    # --- Phase 6: merge causal + extra per chromosome, load genotypes ---
    G_blocks = []
    chr_labels_all = []
    pos_labels_all = []
    causal_mask_all = []

    for chrom in sorted_chroms:
        positions = chr_positions[chrom]

        causal_local_set = causal_by_chr.get(chrom, set())
        extra_local_set = extra_by_chr.get(chrom, set())
        combined_local = np.sort(np.array(list(causal_local_set | extra_local_set),
                                          dtype=int))

        if len(combined_local) == 0:
            continue

        is_causal = np.isin(combined_local, np.array(list(causal_local_set), dtype=int)) \
            if len(causal_local_set) > 0 else np.zeros(len(combined_local), dtype=bool)

        n_combined = len(combined_local)
        n_causal_chr = is_causal.sum()
        print(f"    Chr {chrom:>2}: {n_causal_chr} causal + "
              f"{n_combined - n_causal_chr} non-causal = {n_combined} total")

        # Read genotypes
        with open_bed(chr_files[chrom]) as bed:
            G_chr = bed.read(
                index=(sample_idx.tolist(), combined_local.tolist()), dtype='float64'
            )
        G_chr = np.nan_to_num(G_chr, nan=0.0)

        G_blocks.append(G_chr)
        chr_labels_all.extend([chrom] * n_combined)
        pos_labels_all.extend(positions[combined_local].tolist())
        causal_mask_all.extend(is_causal.tolist())

    # Concatenate
    G_all = np.hstack(G_blocks)
    snp_chrs = np.array(chr_labels_all, dtype=str)
    snp_pos = np.array(pos_labels_all)
    causal_mask = np.array(causal_mask_all, dtype=bool)

    print(f"  Combined: {G_all.shape[0]} samples x {G_all.shape[1]} SNPs "
          f"({causal_mask.sum()} causal)")

    return _finalise(G_all, snp_chrs, snp_pos, causal_mask, sample_idx)


def _finalise(G_all, snp_chrs, snp_pos, causal_mask, sample_idx):
    """Standardise genotypes, build info dataframes."""
    # Standardise each SNP
    G_mean = G_all.mean(axis=0)
    G_std = G_all.std(axis=0)
    G_std[G_std == 0] = 1.0
    G_all = (G_all - G_mean) / G_std

    G_causal = G_all[:, causal_mask]

    snp_info = pd.DataFrame({
        "chr": snp_chrs,
        "pos": snp_pos,
        "is_causal": causal_mask,
    })
    causal_info = snp_info[snp_info["is_causal"]].reset_index(drop=True)

    n_total = G_all.shape[1]
    n_causal = causal_mask.sum()
    unique_chrs = sorted(set(snp_chrs), key=lambda c: int(c) if c.isdigit() else 99)
    print(f"  Final: {n_total} PRS SNPs, {n_causal} causal, "
          f"{n_total - n_causal} non-causal")
    print(f"  Chromosomes: {unique_chrs}")
    print(f"  Ratio: {n_total / n_causal:.0f} SNPs per causal variant")

    return G_all, G_causal, causal_mask, snp_info, causal_info, sample_idx


# ---------------------------------------------------------------------------
# 2. Phenotype simulation — uses ONLY causal SNPs
# ---------------------------------------------------------------------------

def simulate_phenotype(G_causal, gamma, h2=0.2, seed=None):
    """
    Two-pathway interaction model using causal genotypes only.

    Returns: y, P1, P2, set1_idx, set2_idx (indices into G_causal columns)
    """
    rng = np.random.default_rng(seed)
    N, M = G_causal.shape

    perm = rng.permutation(M)
    set1 = perm[:M // 2]
    set2 = perm[M // 2:]

    beta1 = rng.normal(0, 1, len(set1))
    beta2 = rng.normal(0, 1, len(set2))

    g1 = G_causal[:, set1] @ beta1
    g2 = G_causal[:, set2] @ beta2

    # Scale to target h2/2 per pathway
    target_var = h2 / 2.0
    g1 = g1 / g1.std() * np.sqrt(target_var)
    g2 = g2 / g2.std() * np.sqrt(target_var)

    noise_var = 1.0 - target_var
    P1 = g1 + rng.normal(0, np.sqrt(noise_var), N)
    P2 = g2 + rng.normal(0, np.sqrt(noise_var), N)

    P1 = (P1 - P1.mean()) / P1.std()
    P2 = (P2 - P2.mean()) / P2.std()

    y = P1 + P2 + gamma * P1 * P2
    y = (y - y.mean()) / y.std()

    return y, P1, P2, set1, set2


# ---------------------------------------------------------------------------
# 3. GWAS on ALL SNPs
# ---------------------------------------------------------------------------

def compute_gwas(G_train, y_train):
    """Marginal OLS GWAS on every column of G_train."""
    N, P = G_train.shape

    # Vectorised correlation-based GWAS
    y_z = (y_train - y_train.mean()) / y_train.std()
    G_z = G_train  # already standardised

    r = (G_z.T @ y_z) / N
    r2 = r ** 2
    r2 = np.clip(r2, 0, 1 - 1e-10)
    t_stat = r * np.sqrt((N - 2) / (1.0 - r2))
    pvals = 2 * stats.t.sf(np.abs(t_stat), N - 2)
    betas = r

    return betas, pvals


# ---------------------------------------------------------------------------
# 4. PRS on ALL SNPs (with optional p-value threshold)
# ---------------------------------------------------------------------------

def build_prs(G_test, betas, pvals=None, p_threshold=None):
    """
    Build PRS = G_test @ betas, using only SNPs with p < p_threshold if set.
    Standardises the result.
    """
    if p_threshold is not None and pvals is not None:
        mask = pvals < p_threshold
        if mask.sum() == 0:
            return np.random.default_rng(0).standard_normal(G_test.shape[0])
        prs = G_test[:, mask] @ betas[mask]
    else:
        prs = G_test @ betas

    if prs.std() > 0:
        prs = (prs - prs.mean()) / prs.std()
    return prs


def build_chr_prs(G_test, betas, snp_info, pvals=None, p_threshold=None):
    """Build per-chromosome PRS on ALL SNPs for pc-CE."""
    chrs = snp_info["chr"].values
    chr_prs = {}

    for c in sorted(set(chrs)):
        chr_mask = chrs == c
        if p_threshold is not None and pvals is not None:
            sig_mask = pvals < p_threshold
            combined = chr_mask & sig_mask
        else:
            combined = chr_mask

        if combined.sum() == 0:
            continue

        prs_c = G_test[:, combined] @ betas[combined]
        if prs_c.std() > 0:
            prs_c = (prs_c - prs_c.mean()) / prs_c.std()
            chr_prs[c] = prs_c

    return chr_prs


# ---------------------------------------------------------------------------
# 5. CE tests
# ---------------------------------------------------------------------------

def wg_ce_test(y, prs):
    """wg-CE: y ~ alpha * PRS + gamma * PRS^2. Returns gamma_hat, se, p."""
    import statsmodels.api as sm
    X = sm.add_constant(np.column_stack([prs, prs ** 2]))
    model = sm.OLS(y, X).fit()
    return model.params[2], model.bse[2], model.pvalues[2]


def pc_ce_test(y, chr_prs):
    """pc-CE: joint F-test on all chr_i x chr_j interactions."""
    import statsmodels.api as sm

    keys = sorted(chr_prs.keys())
    if len(keys) < 2:
        return np.nan, np.nan, 1.0

    add_cols = [chr_prs[c] for c in keys]
    int_cols = []
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            int_cols.append(chr_prs[keys[i]] * chr_prs[keys[j]])

    if len(int_cols) == 0:
        return np.nan, np.nan, 1.0

    X_add = np.column_stack(add_cols)
    X_int = np.column_stack(int_cols)
    X_full = sm.add_constant(np.column_stack([X_add, X_int]))
    X_red = sm.add_constant(X_add)

    m_full = sm.OLS(y, X_full).fit()
    m_red = sm.OLS(y, X_red).fit()

    n_int = X_int.shape[1]
    N = len(y)
    df_full = N - X_full.shape[1]
    f_stat = ((m_red.ssr - m_full.ssr) / n_int) / (m_full.ssr / df_full)
    f_pval = 1 - stats.f.cdf(f_stat, n_int, df_full)

    gammas = m_full.params[1 + X_add.shape[1]:]
    return gammas.mean(), gammas.std() / np.sqrt(len(gammas)), f_pval


# ---------------------------------------------------------------------------
# 6. Main simulation loop
# ---------------------------------------------------------------------------

def run_simulation(G_all, G_causal, causal_mask, snp_info, causal_info,
                   n_gwas, n_test,
                   gamma_values, n_sim=100, h2=0.2, p_threshold=None, seed=42):
    """
    For each (gamma, replicate):
      1. Simulate phenotype from causal SNPs (all N_gwas + N_test individuals)
      2. GWAS on ALL SNPs using the first N_gwas individuals
      3. PRS from ALL SNPs in the independent test sample (next N_test individuals)
      4. wg-CE and pc-CE on the test sample
      5. Oracle CE on true liabilities in the test sample
    """
    results = []

    for gi, true_gamma in enumerate(gamma_values):
        print(f"\n  Gamma = {true_gamma:+.2f} ({gi + 1}/{len(gamma_values)})")

        for sim_i in range(n_sim):
            sim_seed = seed + gi * 10000 + sim_i

            # Simulate phenotype using CAUSAL SNPs only
            y, P1, P2, s1, s2 = simulate_phenotype(
                G_causal, gamma=true_gamma, h2=h2, seed=sim_seed
            )

            # Non-overlapping GWAS and test samples
            y_gwas, y_test = y[:n_gwas], y[n_gwas:n_gwas + n_test]
            G_all_gwas = G_all[:n_gwas]
            G_all_test = G_all[n_gwas:n_gwas + n_test]

            # GWAS on ALL SNPs in GWAS sample
            betas, pvals = compute_gwas(G_all_gwas, y_gwas)

            # Number of nominally significant SNPs
            n_sig = (pvals < 0.05).sum()
            n_causal_sig = (pvals[causal_mask] < 0.05).sum()

            # --- wg-CE on all-SNP PRS in test sample ---
            prs_wg = build_prs(G_all_test, betas, pvals, p_threshold)
            total_liability_test = P1[n_gwas:n_gwas + n_test] + P2[n_gwas:n_gwas + n_test]
            prs_r2 = np.corrcoef(prs_wg, total_liability_test)[0, 1] ** 2
            wg_gamma, wg_se, wg_pval = wg_ce_test(y_test, prs_wg)

            # --- pc-CE on all-SNP per-chr PRS in test sample ---
            chr_prs = build_chr_prs(G_all_test, betas, snp_info, pvals, p_threshold)
            pc_gamma, pc_se, pc_pval = pc_ce_test(y_test, chr_prs)

            # --- Oracle CE (true liabilities in test sample) ---
            P1t = P1[n_gwas:n_gwas + n_test]
            P2t = P2[n_gwas:n_gwas + n_test]
            P1t = (P1t - P1t.mean()) / P1t.std()
            P2t = (P2t - P2t.mean()) / P2t.std()
            total_t = P1t + P2t
            total_t = (total_t - total_t.mean()) / total_t.std()
            oracle_gamma, oracle_se, oracle_pval = wg_ce_test(y_test, total_t)

            results.append({
                "true_gamma": true_gamma,
                "sim": sim_i,
                "prs_r2": prs_r2,
                "n_sig_snps": n_sig,
                "n_causal_sig": n_causal_sig,
                "wg_gamma_hat": wg_gamma,
                "wg_se": wg_se,
                "wg_pval": wg_pval,
                "pc_gamma_hat": pc_gamma,
                "pc_se": pc_se,
                "pc_pval": pc_pval,
                "oracle_gamma_hat": oracle_gamma,
                "oracle_se": oracle_se,
                "oracle_pval": oracle_pval,
            })

            if (sim_i + 1) % 25 == 0:
                print(f"    {sim_i + 1}/{n_sim} done  "
                      f"[sig SNPs: {n_sig}, PRS R2: {prs_r2:.4f}]")

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# 7. Plotting (3-panel layout)
# ---------------------------------------------------------------------------

def plot_results(df, out_prefix, snp_info=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    gamma_vals = sorted(df["true_gamma"].unique())

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    w = 0.2

    # --- Panel A: Oracle gamma recovery ---
    ax = axes[0]
    oracle_med = [df[df["true_gamma"] == g]["oracle_gamma_hat"].median()
                  for g in gamma_vals]
    oracle_q25 = [df[df["true_gamma"] == g]["oracle_gamma_hat"].quantile(0.25)
                  for g in gamma_vals]
    oracle_q75 = [df[df["true_gamma"] == g]["oracle_gamma_hat"].quantile(0.75)
                  for g in gamma_vals]

    ax.errorbar(np.arange(len(gamma_vals)), oracle_med,
                yerr=[np.array(oracle_med) - np.array(oracle_q25),
                      np.array(oracle_q75) - np.array(oracle_med)],
                fmt="D", color="forestgreen", capsize=4, markersize=6,
                linewidth=1.3, label="Oracle (true liability)")
    ax.plot(np.arange(len(gamma_vals)), gamma_vals,
            color="gray", linewidth=1, linestyle="--", alpha=0.5,
            label="Identity (perfect recovery)")
    ax.set_xticks(np.arange(len(gamma_vals)))
    ax.set_xticklabels([f"{g:.1f}" for g in gamma_vals], fontsize=9)
    ax.set_xlabel("Simulated gamma", fontsize=11)
    ax.set_ylabel("Oracle gamma-hat", fontsize=11)
    ax.set_title("A. Oracle gamma recovery", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, loc="upper left")
    ax.axhline(0, color="gray", linewidth=0.5)

    # --- Panel B: -log10(P) for wg-CE and pc-CE ---
    ax = axes[1]
    pos_l = np.arange(len(gamma_vals)) - w / 2
    pos_r = np.arange(len(gamma_vals)) + w / 2

    for method, pos, color, label in [
        ("wg_pval", pos_l, "steelblue", "wg-CE"),
        ("pc_pval", pos_r, "indianred", "pc-CE"),
    ]:
        logp = [df[df["true_gamma"] == g][method].apply(
            lambda x: -np.log10(max(x, 1e-300))).median()
            for g in gamma_vals]
        ax.bar(pos, logp, width=w * 0.9, color=color, alpha=0.8, label=label)

    ax.set_xticks(np.arange(len(gamma_vals)))
    ax.set_xticklabels([f"{g:.1f}" for g in gamma_vals], fontsize=9)
    ax.set_xlabel("Simulated gamma", fontsize=11)
    ax.set_ylabel("-log10(P value)", fontsize=11)
    ax.set_title("B. Significance (clumped PRS)", fontsize=12, fontweight="bold")
    ax.axhline(-np.log10(0.05), color="gray", linewidth=1, linestyle="--",
               label="P = 0.05")
    ax.legend(fontsize=9)

    # --- Panel C: wg-CE vs pc-CE gamma-hat (paired bars) ---
    ax = axes[2]
    non_null = [g for g in gamma_vals if g != 0]
    x_pos = np.arange(len(non_null))

    wg_m = [df[df["true_gamma"] == g]["wg_gamma_hat"].mean() for g in non_null]
    pc_m = [df[df["true_gamma"] == g]["pc_gamma_hat"].mean() for g in non_null]
    n_rep = df[df["true_gamma"] == non_null[0]].shape[0] if non_null else 1
    wg_se = [df[df["true_gamma"] == g]["wg_gamma_hat"].std() / np.sqrt(n_rep)
             for g in non_null]
    pc_se = [df[df["true_gamma"] == g]["pc_gamma_hat"].std() / np.sqrt(n_rep)
             for g in non_null]

    bw = 0.35
    ax.bar(x_pos - bw / 2, wg_m, bw, yerr=wg_se,
           color="steelblue", alpha=0.8, capsize=3, label="wg-CE")
    ax.bar(x_pos + bw / 2, pc_m, bw, yerr=pc_se,
           color="indianred", alpha=0.8, capsize=3, label="pc-CE")
    ax.set_xticks(x_pos)
    ax.set_xticklabels([f"{g:.1f}" for g in non_null], fontsize=9)
    ax.set_xlabel("Simulated gamma", fontsize=11)
    ax.set_ylabel("Mean gamma-hat (clumped PRS)", fontsize=11)
    ax.set_title("C. wg-CE vs pc-CE gamma-hat", fontsize=12, fontweight="bold")
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.legend(fontsize=9)

    for i, g in enumerate(non_null):
        if abs(wg_m[i]) > 1e-6 and abs(pc_m[i]) > 1e-6:
            ratio = abs(wg_m[i] / pc_m[i])
            y_a = max(abs(wg_m[i]), abs(pc_m[i])) * (1.3 if wg_m[i] > 0 else -1.3)
            ax.annotate(f"{ratio:.0f}x", (x_pos[i], y_a),
                        ha="center", fontsize=7, color="gray")

    plt.tight_layout()
    plt.savefig(f"{out_prefix}_gamma_recovery.pdf", dpi=300, bbox_inches="tight")
    plt.savefig(f"{out_prefix}_gamma_recovery.png", dpi=300, bbox_inches="tight")
    print(f"  Saved: {out_prefix}_gamma_recovery.pdf/png")
    plt.close()

    # --- Summary ---
    mean_r2 = df["prs_r2"].mean()
    mean_sig = df["n_sig_snps"].mean()
    mean_causal_sig = df["n_causal_sig"].mean()

    n_total_snps = len(snp_info) if snp_info is not None else "?"
    n_causal = snp_info["is_causal"].sum() if snp_info is not None else "?"

    print(f"\n{'='*80}")
    print(f"SUMMARY  (clumped-PRS design: {n_total_snps} PRS SNPs, "
          f"{n_causal} causal)")
    print(f"{'='*80}")
    print(f"  Mean PRS R2:                    {mean_r2:.4f}")
    print(f"  Mean nominally sig SNPs (p<.05): {mean_sig:.0f}")
    print(f"  Mean causal SNPs among sig:      {mean_causal_sig:.0f}")

    print(f"\nType I error (gamma=0, P < 0.05):")
    null = df[df["true_gamma"] == 0.0]
    if len(null) > 0:
        print(f"  Oracle:  {(null['oracle_pval'] < 0.05).mean():.3f}")
        print(f"  wg-CE:   {(null['wg_pval'] < 0.05).mean():.3f}")
        print(f"  pc-CE:   {(null['pc_pval'] < 0.05).mean():.3f}")

    print(f"\n{'gamma':>7} | {'Oracle g':>10} {'wg-CE g':>10} {'pc-CE g':>10}"
          f" | {'wg/oracle':>10} {'pc/oracle':>10}"
          f" | {'wg bias':>10} {'pc bias':>10}")
    print("-" * 100)
    for g in gamma_vals:
        sub = df[df["true_gamma"] == g]
        wg = sub["wg_gamma_hat"].mean()
        pc = sub["pc_gamma_hat"].mean()
        orc = sub["oracle_gamma_hat"].mean()
        if g != 0 and abs(orc) > 1e-6:
            print(f"  {g:+.1f}  | {orc:>+10.4f} {wg:>+10.4f} {pc:>+10.4f}"
                  f" | {wg/orc:>10.3f} {pc/orc:>10.3f}"
                  f" | {wg-g:>+10.4f} {pc-g:>+10.4f}")
        else:
            print(f"  {g:+.1f}  | {orc:>+10.4f} {wg:>+10.4f} {pc:>+10.4f}"
                  f" | {'N/A':>10} {'N/A':>10}"
                  f" | {wg-g:>+10.4f} {pc-g:>+10.4f}")


# ---------------------------------------------------------------------------
# 8. Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="wg-CE LD simulation v2 — clumped-PRS design"
    )
    bed_group = parser.add_mutually_exclusive_group(required=True)
    bed_group.add_argument("--bed", default=None,
                           help="Path to single plink .bed file")
    bed_group.add_argument("--bed-prefix", default=None,
                           help="Prefix for per-chromosome .bed files "
                                "(e.g. /path/to/ukb_chr -> expects ukb_chr1.bed, ...)")
    parser.add_argument("--chromosomes", type=str, default=None,
                        help="Comma-separated chromosome numbers to load "
                             "(default: 1-22). Only used with --bed-prefix.")
    parser.add_argument("--out", default="results_LD_clumped",
                        help="Output prefix")
    parser.add_argument("--n-gwas", type=int, default=10000,
                        help="Number of individuals for GWAS (default: 10000)")
    parser.add_argument("--n-test", type=int, default=10000,
                        help="Number of non-overlapping individuals for PRS/CE testing "
                             "(default: 10000)")
    parser.add_argument("--n-causal", type=int, default=1000,
                        help="Number of causal SNPs (default: 1000)")
    parser.add_argument("--n-prs-snps", type=int, default=100000,
                        help="Total number of SNPs for PRS (including causal; "
                             "default: 100000)")
    parser.add_argument("--h2", type=float, default=0.2,
                        help="Additive heritability (default: 0.2)")
    parser.add_argument("--n-sim", type=int, default=100,
                        help="Replicates per gamma (default: 100)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed (default: 42)")
    parser.add_argument("--p-threshold", type=float, default=None,
                        help="P-value threshold for PRS SNP inclusion (default: use all)")
    parser.add_argument("--gamma-values", type=str,
                        default="-1.0,-0.5,-0.3,-0.1,0.0,0.1,0.3,0.5,1.0",
                        help="Comma-separated gamma values")

    args = parser.parse_args()
    gamma_values = [float(x) for x in args.gamma_values.split(",")]

    chromosomes = None
    if args.chromosomes:
        chromosomes = [c.strip() for c in args.chromosomes.split(",")]

    n_total = args.n_gwas + args.n_test

    print("=" * 60)
    print("wg-CE simulation v2 — clumped-PRS design")
    print("=" * 60)
    if args.bed:
        print(f"  BED file:        {args.bed}")
    else:
        print(f"  BED prefix:      {args.bed_prefix}")
        print(f"  Chromosomes:     {chromosomes or '1-22'}")
    print(f"  N GWAS:          {args.n_gwas}")
    print(f"  N test:          {args.n_test}")
    print(f"  N total:         {n_total}")
    print(f"  Causal SNPs:     {args.n_causal}")
    print(f"  PRS SNPs:        {args.n_prs_snps}")
    print(f"  h2:              {args.h2}")
    print(f"  N simulations:   {args.n_sim}")
    print(f"  P threshold:     {args.p_threshold or 'all SNPs'}")
    print(f"  Gamma values:    {gamma_values}")
    print(f"  Seed:            {args.seed}")
    print(f"  Output prefix:   {args.out}")
    print()

    t0 = time.time()

    # Load genotypes (n_gwas + n_test individuals)
    print("Step 1: Loading genotypes...")
    G_all, G_causal, causal_mask, snp_info, causal_info, sample_idx = \
        load_genotypes(
            bed_path=args.bed, bed_prefix=args.bed_prefix,
            chromosomes=chromosomes,
            n_samples=n_total,
            n_causal=args.n_causal,
            n_prs_snps=args.n_prs_snps,
            maf_min=0.05,
            seed=args.seed
        )

    # Run simulations
    print("\nStep 2: Running simulations...")
    df = run_simulation(
        G_all, G_causal, causal_mask, snp_info, causal_info,
        n_gwas=args.n_gwas, n_test=args.n_test,
        gamma_values=gamma_values, n_sim=args.n_sim, h2=args.h2,
        p_threshold=args.p_threshold, seed=args.seed
    )

    # Save
    csv_path = f"{args.out}.csv"
    df.to_csv(csv_path, index=False)
    print(f"\n  Results saved: {csv_path}")

    # Plot
    print("\nStep 3: Generating plots...")
    plot_results(df, args.out, snp_info=snp_info)

    elapsed = time.time() - t0
    print(f"\nDone! Total time: {elapsed / 60:.1f} minutes")


if __name__ == "__main__":
    main()
