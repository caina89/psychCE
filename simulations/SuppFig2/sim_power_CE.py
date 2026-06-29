#!/usr/bin/env python3
"""
Power simulation — CE power using clumped-PRS design with realistic LD.

Addresses Reviewer 2 Major Comment 3 with realistic LD from real genotypes

Design:
  - Load per-chromosome non-LD-pruned .bed files (real LD structure preserved)
  - Select M causal SNPs randomly across all chromosomes
    (causal SNPs may be in LD with each other — this is realistic)
  - Select P total SNPs (default 100,000) including all M causal SNPs
    (simulates an LD-clumped PRS: mostly independent, but not perfectly)
  - Load 2 * max(sample_sizes) individuals for non-overlapping GWAS/test design
  - For each (N, gamma) combination:
      * Subsample 2N individuals, simulate phenotype from causal SNPs
      * Use first N for GWAS on ALL P SNPs
      * Use second N (non-overlapping) for PRS construction and CE testing
      * Run wg-CE, pc-CE, cd-CE
      * Also build "PA-FGRS proxy" scores at target R2 values for comparison
  - The PRS R2 is a MEASURED outcome (determined by N_gwas, LD, h2),
    while the PA-FGRS proxy R2 is a controlled parameter.

This gives two complementary views:
  (A) Power with realistic PRS (varies with N through GWAS power)
  (B) Power with PA-FGRS-like scores (varies with target R2)

Usage:
    python sim_power_CE.py \
        --bed-prefix '/path/to/ukb_imp_v3_chr{CHR}.qced' \
        --chromosomes 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
        --out power_sim_v2 \
        --n-causal 1000 --n-prs-snps 100000 \
        --sample-sizes 5000,10000,20000,50000 \
        --n-sim 200 --seed 42

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
# 1. Genotype loading — clumped-PRS design (no LD-window logic)
# ---------------------------------------------------------------------------

def _find_chr_bed(prefix, chrom):
    """
    Find the .bed file for a given chromosome.

    If prefix contains '{CHR}', replace it with the chromosome number.
    Otherwise try common suffix patterns.
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
                   n_samples=100000, n_causal=1000, n_prs_snps=100000,
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

        # Pick causal SNPs
        causal_in_common = np.sort(rng.choice(n_common, size=n_causal, replace=False))
        causal_global = common_indices[causal_in_common]

        # Pick remaining PRS SNPs
        n_extra = n_prs_snps - n_causal
        non_causal_common = np.setdiff1d(np.arange(n_common), causal_in_common)
        extra_in_common = np.sort(rng.choice(len(non_causal_common), size=n_extra, replace=False))
        extra_global = common_indices[non_causal_common[extra_in_common]]

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

    return _finalise(G_all, snp_chrs, snp_pos, causal_mask)


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
    chr_common_local = {}
    chr_positions = {}

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

    # --- Phase 4: pick causal SNPs (genome-wide) ---
    pool_chroms = []
    pool_local_idx = []
    for chrom in sorted_chroms:
        for idx in chr_common_local[chrom]:
            pool_chroms.append(chrom)
            pool_local_idx.append(idx)

    pool_chroms = np.array(pool_chroms)
    pool_local_idx = np.array(pool_local_idx)

    causal_pick = rng.choice(total_common, size=n_causal, replace=False)
    causal_chroms = pool_chroms[causal_pick]
    causal_local_indices = pool_local_idx[causal_pick]

    causal_by_chr = {}
    for i, chrom in enumerate(causal_chroms):
        causal_by_chr.setdefault(chrom, set()).add(causal_local_indices[i])

    n_causal_per_chr = {c: len(v) for c, v in sorted(causal_by_chr.items(),
                        key=lambda x: int(x[0]) if x[0].isdigit() else 99)}
    print(f"  Causal SNPs per chromosome: {n_causal_per_chr}")
    print(f"  Total causal: {sum(n_causal_per_chr.values())}")

    # --- Phase 5: pick remaining PRS SNPs ---
    n_extra = n_prs_snps - n_causal
    non_causal_pool = np.setdiff1d(np.arange(total_common), causal_pick)
    extra_pick = rng.choice(len(non_causal_pool), size=n_extra, replace=False)
    extra_pool_idx = non_causal_pool[extra_pick]
    extra_chroms = pool_chroms[extra_pool_idx]
    extra_local_indices = pool_local_idx[extra_pool_idx]

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

        with open_bed(chr_files[chrom]) as bed:
            G_chr = bed.read(
                index=(sample_idx.tolist(), combined_local.tolist()), dtype='float64'
            )
        G_chr = np.nan_to_num(G_chr, nan=0.0)

        G_blocks.append(G_chr)
        chr_labels_all.extend([chrom] * n_combined)
        pos_labels_all.extend(positions[combined_local].tolist())
        causal_mask_all.extend(is_causal.tolist())

    G_all = np.hstack(G_blocks)
    snp_chrs = np.array(chr_labels_all, dtype=str)
    snp_pos = np.array(pos_labels_all)
    causal_mask = np.array(causal_mask_all, dtype=bool)

    print(f"  Combined: {G_all.shape[0]} samples x {G_all.shape[1]} SNPs "
          f"({causal_mask.sum()} causal)")

    return _finalise(G_all, snp_chrs, snp_pos, causal_mask)


def _finalise(G_all, snp_chrs, snp_pos, causal_mask):
    """Standardise genotypes, build info dataframes."""
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

    n_total = G_all.shape[1]
    n_causal = causal_mask.sum()
    unique_chrs = sorted(set(snp_chrs), key=lambda c: int(c) if c.isdigit() else 99)
    print(f"  Final: {n_total} PRS SNPs, {n_causal} causal, "
          f"{n_total - n_causal} non-causal")
    print(f"  Chromosomes: {unique_chrs}")
    print(f"  Ratio: {n_total / n_causal:.0f} SNPs per causal variant")

    return G_all, G_causal, causal_mask, snp_info


# ---------------------------------------------------------------------------
# 2. Phenotype simulation (two-pathway model, causal SNPs only)
# ---------------------------------------------------------------------------

def simulate_phenotype(G_causal, gamma, h2=0.2, seed=None):
    """
    Two-pathway interaction model. Returns y, P1, P2.
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

    return y, P1, P2


# ---------------------------------------------------------------------------
# 3. GWAS (vectorised, on all SNPs)
# ---------------------------------------------------------------------------

def compute_gwas(G_train, y_train):
    """Marginal OLS GWAS on every column of G_train. Fast vectorised version."""
    N = G_train.shape[0]
    y_z = (y_train - y_train.mean()) / y_train.std()

    r = (G_train.T @ y_z) / N
    r2 = np.clip(r ** 2, 0, 1 - 1e-10)
    t_stat = r * np.sqrt((N - 2) / (1.0 - r2))
    pvals = 2 * stats.t.sf(np.abs(t_stat), N - 2)

    return r, pvals


# ---------------------------------------------------------------------------
# 4. PRS construction
# ---------------------------------------------------------------------------

def build_prs_allsnp(G_test, betas, pvals=None, p_threshold=None):
    """PRS from ALL SNPs (or p-value filtered). Returns standardised score."""
    if p_threshold is not None and pvals is not None:
        mask = pvals < p_threshold
        if mask.sum() == 0:
            return np.zeros(G_test.shape[0])
        prs = G_test[:, mask] @ betas[mask]
    else:
        prs = G_test @ betas

    if prs.std() > 0:
        prs = (prs - prs.mean()) / prs.std()
    return prs


def build_chr_prs_allsnp(G_test, betas, snp_info, pvals=None, p_threshold=None):
    """Per-chromosome PRS from ALL SNPs."""
    chrs = snp_info["chr"].values
    chr_prs = {}
    for c in sorted(set(chrs)):
        cmask = chrs == c
        if p_threshold is not None and pvals is not None:
            cmask = cmask & (pvals < p_threshold)
        if cmask.sum() == 0:
            continue
        s = G_test[:, cmask] @ betas[cmask]
        if s.std() > 0:
            s = (s - s.mean()) / s.std()
            chr_prs[c] = s
    return chr_prs


def build_score_at_r2(true_liability, target_r2, N, rng):
    """
    Artificial score at controlled R2 (models PA-FGRS).
    score = sqrt(R2)*Z_true + sqrt(1-R2)*Z_noise
    """
    z_true = (true_liability - true_liability.mean()) / true_liability.std()
    z_noise = rng.standard_normal(N)
    z_noise = (z_noise - z_noise.mean()) / z_noise.std()
    score = np.sqrt(target_r2) * z_true + np.sqrt(1.0 - target_r2) * z_noise
    score = (score - score.mean()) / score.std()
    return score


# ---------------------------------------------------------------------------
# 5. CE tests
# ---------------------------------------------------------------------------

def wg_ce_test(y, score):
    """wg-CE: y ~ a*S + g*S^2. Returns gamma_hat, p_value."""
    import statsmodels.api as sm
    X = sm.add_constant(np.column_stack([score, score ** 2]))
    m = sm.OLS(y, X).fit()
    return m.params[2], m.pvalues[2]


def pc_ce_test(y, scores_dict):
    """pc-CE: joint F-test on chr_i x chr_j interactions."""
    import statsmodels.api as sm

    keys = sorted(scores_dict.keys())
    if len(keys) < 2:
        return np.nan, 1.0

    add_cols = [scores_dict[k] for k in keys]
    int_cols = []
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            int_cols.append(scores_dict[keys[i]] * scores_dict[keys[j]])
    if len(int_cols) == 0:
        return np.nan, 1.0

    X_add = np.column_stack(add_cols)
    X_int = np.column_stack(int_cols)
    X_full = sm.add_constant(np.column_stack([X_add, X_int]))
    X_red = sm.add_constant(X_add)

    m_full = sm.OLS(y, X_full).fit()
    m_red = sm.OLS(y, X_red).fit()

    n_int = X_int.shape[1]
    df_full = len(y) - X_full.shape[1]
    f_stat = ((m_red.ssr - m_full.ssr) / n_int) / (m_full.ssr / df_full)
    f_pval = 1 - stats.f.cdf(f_stat, n_int, df_full)

    gammas = m_full.params[1 + X_add.shape[1]:]
    return gammas.mean(), f_pval


def cd_ce_test(y, score1, score2):
    """cd-CE: y ~ a1*S1 + a2*S2 + g*S1*S2. Returns gamma_hat, p_value."""
    import statsmodels.api as sm
    X = sm.add_constant(np.column_stack([score1, score2, score1 * score2]))
    m = sm.OLS(y, X).fit()
    return m.params[3], m.pvalues[3]


# ---------------------------------------------------------------------------
# 6. Main simulation
# ---------------------------------------------------------------------------

def run_simulation(G_all, G_causal, causal_mask, snp_info,
                   sample_sizes, gamma_values, fgrs_r2_values,
                   n_sim=200, h2=0.2, p_threshold=None, seed=42):
    """
    Two-part power analysis with non-overlapping GWAS and test samples.

    Part A -- Real PRS (with LD):
      For each (N, gamma): subsample 2N individuals, simulate phenotype for
      all 2N, use first N for GWAS, second N for PRS construction and CE tests.
      The PRS R2 is a measured outcome.

    Part B -- Modelled PA-FGRS:
      For each (N, gamma, target_R2): subsample 2N individuals, simulate
      phenotype for all 2N, use second N for CE tests with artificial scores
      at target R2 on the true liability.
      Uses the largest available N to isolate the effect of score accuracy.

    Returns: (df_prs, df_fgrs) -- two DataFrames.
    """
    N_max = G_all.shape[0]
    rng_master = np.random.default_rng(seed)
    results_prs = []
    results_fgrs = []

    # ===================================================================
    # Part A: Real PRS with LD (N GWAS + N test, non-overlapping)
    # ===================================================================
    total_a = len(sample_sizes) * len(gamma_values)
    combo = 0

    for N in sample_sizes:
        N_total_needed = 2 * N  # N for GWAS + N for test
        if N_total_needed > N_max:
            print(f"  SKIP N={N}: need 2N={N_total_needed} but only {N_max} loaded")
            continue

        for gamma in gamma_values:
            combo += 1
            print(f"  PRS [{combo}/{total_a}] N_gwas={N:,}, N_test={N:,}, gamma={gamma:+.2f}")

            for sim_i in range(n_sim):
                sim_seed = seed + combo * 100000 + sim_i
                rng = np.random.default_rng(sim_seed)

                # Subsample 2N individuals (N GWAS + N test)
                idx = rng.choice(N_max, size=N_total_needed, replace=False)
                G_sub = G_all[idx]
                Gc_sub = G_causal[idx]

                # Simulate phenotype for all 2N individuals
                y, P1, P2 = simulate_phenotype(Gc_sub, gamma, h2, sim_seed)

                # Non-overlapping split
                y_gwas, y_test = y[:N], y[N:]
                G_gwas = G_sub[:N]
                G_test = G_sub[N:]

                # GWAS on ALL SNPs in GWAS sample
                betas, pvals = compute_gwas(G_gwas, y_gwas)

                # wg-CE: whole-genome PRS in test sample
                prs_wg = build_prs_allsnp(G_test, betas, pvals, p_threshold)
                total_liab = P1[N:] + P2[N:]
                prs_r2 = np.corrcoef(prs_wg, total_liab)[0, 1] ** 2

                g_wg, p_wg = wg_ce_test(y_test, prs_wg)

                # pc-CE: per-chromosome PRS in test sample
                chr_prs = build_chr_prs_allsnp(
                    G_test, betas, snp_info, pvals, p_threshold
                )
                g_pc, p_pc = pc_ce_test(y_test, chr_prs)

                # cd-CE: separate GWAS for each pathway, PRS in test sample
                betas_p1, pvals_p1 = compute_gwas(G_gwas, P1[:N])
                betas_p2, pvals_p2 = compute_gwas(G_gwas, P2[:N])
                prs_d1 = build_prs_allsnp(G_test, betas_p1, pvals_p1, p_threshold)
                prs_d2 = build_prs_allsnp(G_test, betas_p2, pvals_p2, p_threshold)
                g_cd, p_cd = cd_ce_test(y_test, prs_d1, prs_d2)

                results_prs.append({
                    "N": N, "true_gamma": gamma, "sim": sim_i,
                    "prs_r2": prs_r2,
                    "wg_power": int(p_wg < 0.05), "wg_gamma": g_wg, "wg_pval": p_wg,
                    "pc_power": int(p_pc < 0.05), "pc_gamma": g_pc, "pc_pval": p_pc,
                    "cd_power": int(p_cd < 0.05), "cd_gamma": g_cd, "cd_pval": p_cd,
                })

            # Progress
            sub = [r for r in results_prs if r["N"] == N and r["true_gamma"] == gamma]
            mean_r2 = np.mean([r["prs_r2"] for r in sub])
            wg_pow = np.mean([r["wg_power"] for r in sub])
            cd_pow = np.mean([r["cd_power"] for r in sub])
            print(f"    R2={mean_r2:.4f}, wg-CE={wg_pow:.1%}, cd-CE={cd_pow:.1%}")

    # ===================================================================
    # Part B: Modelled PA-FGRS (artificial score at target R2)
    # ===================================================================
    # Use the largest N where we can fit 2N individuals
    N_fgrs = max(s for s in sample_sizes if 2 * s <= N_max)
    N_fgrs_total = 2 * N_fgrs  # N for "train" (unused for FGRS) + N for test
    total_b = len(fgrs_r2_values) * len(gamma_values)
    combo = 0

    for target_r2 in fgrs_r2_values:
        for gamma in gamma_values:
            combo += 1
            print(f"  FGRS [{combo}/{total_b}] R2={target_r2:.2f}, gamma={gamma:+.2f}, "
                  f"N_test={N_fgrs:,}")

            for sim_i in range(n_sim):
                sim_seed = seed + 9999 * 100000 + combo * 10000 + sim_i
                rng = np.random.default_rng(sim_seed)

                idx = rng.choice(N_max, size=N_fgrs_total, replace=False)
                Gc_sub = G_causal[idx]

                # Simulate phenotype for all 2N individuals
                y, P1, P2 = simulate_phenotype(Gc_sub, gamma, h2, sim_seed)

                # Use second half as test sample
                y_test = y[N_fgrs:]
                P1_test = P1[N_fgrs:]
                P2_test = P2[N_fgrs:]
                total_liab = P1_test + P2_test

                # wg-CE: single score at target R2
                score_wg = build_score_at_r2(total_liab, target_r2, N_fgrs, rng)
                g_wg, p_wg = wg_ce_test(y_test, score_wg)

                # pc-CE: split score into 10 pseudo-chromosome partitions
                partition_scores = {}
                n_part = 10
                for k in range(n_part):
                    sub_r2 = target_r2 / n_part
                    partition_scores[str(k)] = build_score_at_r2(
                        total_liab, sub_r2, N_fgrs, rng
                    )
                g_pc, p_pc = pc_ce_test(y_test, partition_scores)

                # cd-CE: separate scores for each pathway
                score_d1 = build_score_at_r2(P1_test, target_r2, N_fgrs, rng)
                score_d2 = build_score_at_r2(P2_test, target_r2, N_fgrs, rng)
                g_cd, p_cd = cd_ce_test(y_test, score_d1, score_d2)

                results_fgrs.append({
                    "N": N_fgrs, "R2": target_r2, "true_gamma": gamma, "sim": sim_i,
                    "wg_power": int(p_wg < 0.05), "wg_gamma": g_wg, "wg_pval": p_wg,
                    "pc_power": int(p_pc < 0.05), "pc_gamma": g_pc, "pc_pval": p_pc,
                    "cd_power": int(p_cd < 0.05), "cd_gamma": g_cd, "cd_pval": p_cd,
                })

    df_prs = pd.DataFrame(results_prs)
    df_fgrs = pd.DataFrame(results_fgrs)
    return df_prs, df_fgrs


# ---------------------------------------------------------------------------
# 7. Plotting
# ---------------------------------------------------------------------------

def _format_N(x, pos=None):
    """Format sample size: 5000 -> '5,000'."""
    if x >= 1000:
        return f"{int(x):,}"
    return str(int(x))


def _fix_log_N_axis(ax, n_values):
    """Replace default log ticks with explicit sample-size labels."""
    import matplotlib.ticker as mticker
    ax.set_xticks(sorted(n_values))
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(_format_N))
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.tick_params(axis='x', rotation=45)


def plot_results(df_prs, df_fgrs, out_prefix):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    # Aggregate PRS results
    if len(df_prs) > 0:
        agg_prs = df_prs.groupby(["N", "true_gamma"]).agg(
            wg_power=("wg_power", "mean"),
            pc_power=("pc_power", "mean"),
            cd_power=("cd_power", "mean"),
            prs_r2=("prs_r2", "mean"),
        ).reset_index()

    # Aggregate FGRS results
    if len(df_fgrs) > 0:
        agg_fgrs = df_fgrs.groupby(["R2", "true_gamma"]).agg(
            wg_power=("wg_power", "mean"),
            pc_power=("pc_power", "mean"),
            cd_power=("cd_power", "mean"),
        ).reset_index()

    gamma_nonnull = sorted([g for g in df_prs["true_gamma"].unique() if g != 0]) \
        if len(df_prs) > 0 else []

    # =================================================================
    # FIGURE 1: PRS power vs N, one panel per gamma
    # =================================================================
    if len(df_prs) > 0 and len(gamma_nonnull) > 0:
        n_panels = len(gamma_nonnull)
        fig, axes = plt.subplots(1, n_panels, figsize=(5.5 * n_panels, 5), sharey=True)
        if n_panels == 1:
            axes = [axes]

        for ax, gamma in zip(axes, gamma_nonnull):
            sub = agg_prs[agg_prs["true_gamma"] == gamma].sort_values("N")
            ax.plot(sub["N"], sub["wg_power"], "-o", color="steelblue",
                    markersize=5, linewidth=2, label="wg-CE")
            ax.plot(sub["N"], sub["pc_power"], "--s", color="indianred",
                    markersize=4, linewidth=1.5, alpha=0.7, label="pc-CE")
            ax.plot(sub["N"], sub["cd_power"], ":^", color="forestgreen",
                    markersize=5, linewidth=2, alpha=0.8, label="cd-CE")

            for _, row in sub.iterrows():
                ax.annotate(f"R2={row['prs_r2']:.3f}",
                            (row["N"], max(row["wg_power"], row["cd_power"]) + 0.04),
                            fontsize=7, ha="center", color="gray")

            ax.set_title(f"gamma = {gamma}", fontsize=12)
            ax.set_xlabel("Sample size (N)", fontsize=11)
            ax.axhline(0.8, color="gray", linewidth=0.5, linestyle=":")
            ax.axhline(0.05, color="lightgray", linewidth=0.5, linestyle="--")
            ax.set_ylim(-0.02, 1.05)
            ax.set_xscale("log")
            _fix_log_N_axis(ax, sub["N"].unique())
            ax.legend(fontsize=9)
            ax.tick_params(labelsize=9)

        axes[0].set_ylabel("Power (P < 0.05)", fontsize=11)
        fig.suptitle("Power with clumped PRS (real LD)", fontsize=13, y=1.02)
        plt.tight_layout()
        fig.savefig(f"{out_prefix}_prs_power_vs_N.pdf", dpi=300, bbox_inches="tight")
        fig.savefig(f"{out_prefix}_prs_power_vs_N.png", dpi=300, bbox_inches="tight")
        print(f"  Saved {out_prefix}_prs_power_vs_N.pdf")
        plt.close(fig)

    # =================================================================
    # FIGURE 2: PA-FGRS proxy power vs R2, one panel per gamma
    # =================================================================
    if len(df_fgrs) > 0:
        gamma_nonnull_f = sorted([g for g in df_fgrs["true_gamma"].unique() if g != 0])
        if gamma_nonnull_f:
            n_panels = len(gamma_nonnull_f)
            fig, axes = plt.subplots(1, n_panels, figsize=(5.5 * n_panels, 5), sharey=True)
            if n_panels == 1:
                axes = [axes]

            for ax, gamma in zip(axes, gamma_nonnull_f):
                sub = agg_fgrs[agg_fgrs["true_gamma"] == gamma].sort_values("R2")
                ax.plot(sub["R2"], sub["wg_power"], "-o", color="steelblue",
                        markersize=5, linewidth=2, label="wg-CE")
                ax.plot(sub["R2"], sub["pc_power"], "--s", color="indianred",
                        markersize=4, linewidth=1.5, alpha=0.7, label="pc-CE")
                ax.plot(sub["R2"], sub["cd_power"], ":^", color="forestgreen",
                        markersize=5, linewidth=2, alpha=0.8, label="cd-CE")

                ax.set_title(f"gamma = {gamma}", fontsize=12)
                ax.set_xlabel("Score R2 (PA-FGRS proxy)", fontsize=11)
                ax.axhline(0.8, color="gray", linewidth=0.5, linestyle=":")
                ax.axhline(0.05, color="lightgray", linewidth=0.5, linestyle="--")
                ax.set_ylim(-0.02, 1.05)
                ax.legend(fontsize=9)
                ax.tick_params(labelsize=9)

            axes[0].set_ylabel("Power (P < 0.05)", fontsize=11)
            N_used = int(df_fgrs["N"].iloc[0])
            fig.suptitle(f"Power with PA-FGRS proxy scores (N={N_used:,})",
                         fontsize=13, y=1.02)
            plt.tight_layout()
            fig.savefig(f"{out_prefix}_fgrs_power_vs_R2.pdf", dpi=300, bbox_inches="tight")
            fig.savefig(f"{out_prefix}_fgrs_power_vs_R2.png", dpi=300, bbox_inches="tight")
            print(f"  Saved {out_prefix}_fgrs_power_vs_R2.pdf")
            plt.close(fig)

    # =================================================================
    # FIGURE 3: Type I error (gamma = 0)
    # =================================================================
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    legend_items = [
        Line2D([0], [0], color="red", linestyle="-", linewidth=1, label="Nominal 5%"),
        Line2D([0], [0], color="steelblue", linestyle="-", marker="o", markersize=4, label="wg-CE"),
        Line2D([0], [0], color="indianred", linestyle="--", marker="s", markersize=3, label="pc-CE"),
        Line2D([0], [0], color="forestgreen", linestyle=":", marker="^", markersize=4, label="cd-CE"),
    ]

    # Type I error vs N (PRS)
    ax = axes[0]
    if len(df_prs) > 0:
        null_prs = agg_prs[agg_prs["true_gamma"] == 0].sort_values("N")
        if len(null_prs) > 0:
            ax.plot(null_prs["N"], null_prs["wg_power"], "-o", color="steelblue", markersize=5)
            ax.plot(null_prs["N"], null_prs["pc_power"], "--s", color="indianred", markersize=4, alpha=0.7)
            ax.plot(null_prs["N"], null_prs["cd_power"], ":^", color="forestgreen", markersize=5, alpha=0.8)
    ax.axhline(0.05, color="red", linewidth=1, linestyle="-")
    ax.set_xlabel("Sample size (N)", fontsize=11)
    ax.set_ylabel("Type I error rate", fontsize=11)
    ax.set_title("Type I error vs N (clumped PRS)", fontsize=12)
    ax.set_xscale("log")
    if len(df_prs) > 0:
        _fix_log_N_axis(ax, df_prs["N"].unique())
    ax.set_ylim(-0.01, 0.20)
    ax.legend(handles=legend_items, fontsize=8, loc="upper left")

    # Type I error vs R2 (FGRS)
    ax = axes[1]
    if len(df_fgrs) > 0:
        null_fgrs = agg_fgrs[agg_fgrs["true_gamma"] == 0].sort_values("R2")
        if len(null_fgrs) > 0:
            ax.plot(null_fgrs["R2"], null_fgrs["wg_power"], "-o", color="steelblue", markersize=5)
            ax.plot(null_fgrs["R2"], null_fgrs["pc_power"], "--s", color="indianred", markersize=4, alpha=0.7)
            ax.plot(null_fgrs["R2"], null_fgrs["cd_power"], ":^", color="forestgreen", markersize=5, alpha=0.8)
    ax.axhline(0.05, color="red", linewidth=1, linestyle="-")
    ax.set_xlabel("Score R2 (PA-FGRS proxy)", fontsize=11)
    ax.set_ylabel("Type I error rate", fontsize=11)
    ax.set_title("Type I error vs R2 (PA-FGRS proxy)", fontsize=12)
    ax.set_ylim(-0.01, 0.20)
    ax.legend(handles=legend_items, fontsize=8, loc="upper left")

    plt.tight_layout()
    fig.savefig(f"{out_prefix}_type1error.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out_prefix}_type1error.png", dpi=300, bbox_inches="tight")
    print(f"  Saved {out_prefix}_type1error.pdf")
    plt.close(fig)

    # =================================================================
    # Summary tables
    # =================================================================
    print(f"\n{'='*80}")
    print("PART A: REAL PRS POWER (clumped PRS, N GWAS + N test non-overlapping)")
    print(f"{'='*80}")
    if len(df_prs) > 0:
        print(f"{'N_gwas':>8} {'gamma':>7} | {'R2':>6} {'wg-CE':>8} {'pc-CE':>8} {'cd-CE':>8}")
        print("-" * 60)
        for _, row in agg_prs.sort_values(["true_gamma", "N"]).iterrows():
            print(f"{row['N']:>8,.0f} {row['true_gamma']:>+7.2f} | "
                  f"{row['prs_r2']:>6.4f} "
                  f"{row['wg_power']:>7.1%} {row['pc_power']:>7.1%} {row['cd_power']:>7.1%}")

    print(f"\n{'='*80}")
    print("PART B: PA-FGRS PROXY POWER (modelled score at target R2)")
    print(f"{'='*80}")
    if len(df_fgrs) > 0:
        print(f"{'R2':>6} {'gamma':>7} | {'wg-CE':>8} {'pc-CE':>8} {'cd-CE':>8}")
        print("-" * 50)
        for _, row in agg_fgrs.sort_values(["true_gamma", "R2"]).iterrows():
            print(f"{row['R2']:>6.2f} {row['true_gamma']:>+7.2f} | "
                  f"{row['wg_power']:>7.1%} {row['pc_power']:>7.1%} {row['cd_power']:>7.1%}")


# ---------------------------------------------------------------------------
# 8. Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="CE power simulation v2 — clumped-PRS design"
    )
    bed_group = parser.add_mutually_exclusive_group(required=True)
    bed_group.add_argument("--bed", default=None,
                           help="Path to single plink .bed file")
    bed_group.add_argument("--bed-prefix", default=None,
                           help="Prefix for per-chromosome .bed files "
                                "(e.g. /path/to/ukb_chr -> expects ukb_chr1.bed, ...)")
    parser.add_argument("--chromosomes", type=str, default=None,
                        help="Comma-separated chromosomes (default: 1-22). "
                             "Only with --bed-prefix.")
    parser.add_argument("--out", default="power_sim_v2",
                        help="Output prefix")
    parser.add_argument("--n-causal", type=int, default=1000,
                        help="Causal SNPs (default: 1000)")
    parser.add_argument("--n-prs-snps", type=int, default=100000,
                        help="Total number of SNPs for PRS (including causal; "
                             "default: 100000)")
    parser.add_argument("--h2", type=float, default=0.2,
                        help="Additive heritability (default: 0.2)")
    parser.add_argument("--n-sim", type=int, default=200,
                        help="Replicates per condition (default: 200)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed")
    parser.add_argument("--p-threshold", type=float, default=None,
                        help="P-value threshold for PRS (default: all SNPs)")
    parser.add_argument("--sample-sizes", type=str,
                        default="5000,10000,20000,50000",
                        help="Comma-separated GWAS sample sizes (N). "
                             "For each N, N non-overlapping test individuals are also used "
                             "(total 2N per condition). Default: 5000,10000,20000,50000")
    parser.add_argument("--gamma-values", type=str,
                        default="-0.10,-0.05,-0.03,0.0",
                        help="Comma-separated gamma values")
    parser.add_argument("--fgrs-r2-values", type=str,
                        default="0.01,0.02,0.05,0.10,0.15,0.20",
                        help="Comma-separated R2 values for PA-FGRS proxy")

    args = parser.parse_args()

    sample_sizes = [int(x) for x in args.sample_sizes.split(",")]
    gamma_values = [float(x) for x in args.gamma_values.split(",")]
    fgrs_r2_values = [float(x) for x in args.fgrs_r2_values.split(",")]

    # Need 2 * max(sample_sizes) individuals for non-overlapping GWAS/test
    n_load = 2 * max(sample_sizes)

    chromosomes = None
    if args.chromosomes:
        chromosomes = [c.strip() for c in args.chromosomes.split(",")]

    print("=" * 60)
    print("CE Power Simulation v2 — clumped-PRS design")
    print("=" * 60)
    if args.bed:
        print(f"  BED file:        {args.bed}")
    else:
        print(f"  BED prefix:      {args.bed_prefix}")
        print(f"  Chromosomes:     {chromosomes or '1-22'}")
    print(f"  Causal SNPs:     {args.n_causal}")
    print(f"  PRS SNPs:        {args.n_prs_snps}")
    print(f"  h2:              {args.h2}")
    print(f"  N simulations:   {args.n_sim}")
    print(f"  P threshold:     {args.p_threshold or 'all SNPs'}")
    print(f"  GWAS sample sizes (N): {sample_sizes}")
    print(f"  Test sample sizes:     {sample_sizes}  (non-overlapping)")
    print(f"  Gamma values:    {gamma_values}")
    print(f"  FGRS R2 values:  {fgrs_r2_values}")
    print(f"  Individuals to load: {n_load}  (2 x max N)")
    n_cond_a = len(sample_sizes) * len(gamma_values)
    n_cond_b = len(fgrs_r2_values) * len(gamma_values)
    print(f"  PRS conditions:  {n_cond_a} ({n_cond_a * args.n_sim:,} sims)")
    print(f"  FGRS conditions: {n_cond_b} ({n_cond_b * args.n_sim:,} sims)")
    print()

    t0 = time.time()

    print("Step 1: Loading genotypes...")
    G_all, G_causal, causal_mask, snp_info = load_genotypes(
        bed_path=args.bed, bed_prefix=args.bed_prefix,
        chromosomes=chromosomes,
        n_samples=n_load,
        n_causal=args.n_causal,
        n_prs_snps=args.n_prs_snps,
        maf_min=0.05,
        seed=args.seed
    )

    print("\nStep 2: Running simulations...")
    df_prs, df_fgrs = run_simulation(
        G_all, G_causal, causal_mask, snp_info,
        sample_sizes=sample_sizes,
        gamma_values=gamma_values,
        fgrs_r2_values=fgrs_r2_values,
        n_sim=args.n_sim, h2=args.h2,
        p_threshold=args.p_threshold, seed=args.seed
    )

    # Save
    df_prs.to_csv(f"{args.out}_prs.csv", index=False)
    df_fgrs.to_csv(f"{args.out}_fgrs.csv", index=False)
    print(f"\n  Saved: {args.out}_prs.csv, {args.out}_fgrs.csv")

    print("\nStep 3: Plotting...")
    plot_results(df_prs, df_fgrs, args.out)

    elapsed = time.time() - t0
    print(f"\nDone! Total time: {elapsed / 60:.1f} minutes")


if __name__ == "__main__":
    main()
