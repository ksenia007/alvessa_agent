
"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2025-09-18
Updated: 


Description: 

2 Sample MR tool to perform Mendelian Randomization analysis using genal python package"""

from __future__ import annotations
import os, math, urllib.request
from typing import Any
import pandas as pd
import genal

from src.state import State
from src.tools.base import Node
import os, math, urllib.request
from src.state import State
from src.tools.base import Node

DEBUG = True

# Load manifest once
MANIFEST_PATH = "local_dbs/FinnGen/finngen_R12_manifest.tsv"
AVAIL_GWAS = pd.read_csv(MANIFEST_PATH, sep="\t", low_memory=False)
PHENOCODE_SET = set(AVAIL_GWAS["phenocode"].astype(str))

# Build the “ID: description” catalog for the node description
pairs = AVAIL_GWAS[["phenocode", "phenotype"]].dropna()
id_desc = [f"{row.phenocode}: {row.phenotype}" for row in pairs.itertuples(index=False)]
id_desc_str = "; ".join(id_desc)

def _subset_exposure_to_cis(exp_df: pd.DataFrame,
                            chrom_int: int,
                            start: int,
                            end: int,
                            pad: int = 20_000) -> pd.DataFrame:
    lo = max(0, int(start) - pad)
    hi = int(end) + pad
    return exp_df[(exp_df["#chrom"] == int(chrom_int)) &
                  (exp_df["pos"] >= lo) &
                  (exp_df["pos"] <= hi)]


def build_mr_summary_paragraph(manifest_df: pd.DataFrame,
                               exposure_id: str,
                               outcome_id: str,
                               res: pd.DataFrame,
                               cohort_label: str = "FinnGen R12", 
                               cis_mode: bool = False, 
                               cis_comment: srt = '') -> str:
    # --- manifest meta (guaranteed columns)
    def meta(pid: str):
        row = manifest_df.loc[manifest_df["phenocode"] == pid]
        if row.empty:
            return {"id": pid, "phenotype": "(description unavailable)", "cases": None, "controls": None, "total": None}
        r = row.iloc[0]
        cases = int(r["num_cases"]) if pd.notnull(r["num_cases"]) else None
        ctrls = int(r["num_controls"]) if pd.notnull(r["num_controls"]) else None
        total = (cases or 0) + (ctrls or 0) if (cases is not None or ctrls is not None) else None
        return {"id": pid,
                "phenotype": str(r["phenotype"]) if pd.notnull(r["phenotype"]) else "(description unavailable)",
                "cases": cases, "controls": ctrls, "total": total}

    expm = meta(exposure_id)
    outm = meta(outcome_id)

    # --- pull MR rows
    ivw = res[res["method"].str.contains("Inverse-Variance Weighted")].iloc[0]
    wm  = res[res["method"].str.contains("Weighted Median")].iloc[0]
    eg  = res[res["method"].str.contains("Egger Intercept")].iloc[0]

    # IVW as OR
    b, se, p = ivw["b"], ivw["se"], ivw["pval"]
    OR  = math.exp(b); CIlo = math.exp(b - 1.96*se); CIhi = math.exp(b + 1.96*se)
    direction = "higher" if OR > 1 else "lower"
    signif = "significant" if p < 0.05 else "not significant"

    # Heterogeneity
    q, q_df, q_p = ivw.get("Q", None), ivw.get("Q_df", None), ivw.get("Q_pval", None)
    if q is not None and q_p is not None:
        het_flag = "evidence of heterogeneity" if q_p < 0.05 else "no strong heterogeneity"
        het_txt = f"Heterogeneity: Q={q:.2f}, df={int(q_df)}, p={q_p:.2e} ({het_flag}); heterogeneity means SNP effects vary more than expected by chance and is common with many instruments, so we report IVW random-effects."
    else:
        het_txt = "Heterogeneity: not reported."

    # Weighted Median
    b_wm, se_wm, p_wm = wm["b"], wm["se"], wm["pval"]
    OR_wm, L_wm, H_wm = math.exp(b_wm), math.exp(b_wm - 1.96*se_wm), math.exp(b_wm + 1.96*se_wm)
    wm_consistent = "directionally consistent" if (b * b_wm) > 0 else "in the opposite direction"
    wm_txt = f"Weighted Median was {wm_consistent} (OR={OR_wm:.3f}, 95% CI {L_wm:.3f}–{H_wm:.3f}, p={p_wm:.2e})."

    # MR-Egger intercept
    eg_p = eg["pval"]
    eg_txt = ("MR-Egger intercept not significant (p=" + f"{eg_p:.2e}) suggesting no average directional pleiotropy."
              if eg_p >= 0.05 else
              "MR-Egger intercept significant (p=" + f"{eg_p:.2e}) indicating potential average directional pleiotropy.")

    # Format Ns
    def fmt_n(m):
        parts = []
        if m["total"] is not None: parts.append(f"N={m['total']}")
        if m["cases"] is not None: parts.append(f"cases={m['cases']}")
        if m["controls"] is not None: parts.append(f"controls={m['controls']}")
        return "; " + ", ".join(parts) if parts else ""

    # Final paragraph
    if cis_mode:
        cis_line = cis_comment
    else:
        cis_line = ''
    para = (
        f"Using {cohort_label}, we run 2 sample Mendelian Randomization to test whether genetic liability to {expm['id']} "
        f"({expm['phenotype']}{fmt_n(expm)}) causally influences {outm['id']} "
        f"({outm['phenotype']}{fmt_n(outm)}). {cis_line}"
        "The primary IVW (random-effects) estimate "
        f"indicated {direction} {outcome_id} risk per 1 SD increase in genetically predicted "
        f"{exposure_id}: OR={OR:.3f} (95% CI {CIlo:.3f}–{CIhi:.3f}, p={p:.2e}, {signif}). "
        f"{het_txt} {wm_txt} {eg_txt} Instruments were clumped at p<5e-8, r2=0.1, kb=250 and harmonized with palindromic SNPs handled via allele frequencies (action=2)."
    )
    return para


def _fetch_finngen_gwas(manifest_df: pd.DataFrame, phenocode: str, outdir: str) -> str:
    row = manifest_df.loc[manifest_df["phenocode"] == phenocode]
    if row.empty:
        return None
    url = row["path_https"].values[0]
    os.makedirs(outdir, exist_ok=True)
    local_path = os.path.join(outdir, f"{phenocode}.gz")
    if not os.path.exists(local_path):
        if DEBUG: print(f"[MR] Downloading {url} -> {local_path}")
        urllib.request.urlretrieve(url, local_path)
    else:
        if DEBUG: print(f"[MR] Using cached {local_path}")
    return local_path

def twosample_mr_agent(state: "State",
                       exposure_gwas: str,
                       outcome_gwas: str,
                       gene_name: str | None = None) -> None:
    """
    In-place MR: downloads FinnGen GWAS for exposure/outcome, runs 2-sample MR,
    saves TSV+PNG, and appends a headline summary to state['text_notes'].
    If gene_name is provided and resolvable in state['gene_entities'], runs a cis-only
    analysis using exposure variants within ±20 kb of the gene (hg38).
    Returns None (mutates state only).
    """
    state.setdefault("text_notes", [])
    gene_objs = state.get("gene_entities", {}).copy()

    # Validate phenocodes
    missing = [c for c in (exposure_gwas, outcome_gwas) if c not in PHENOCODE_SET]
    if missing:
        msg = f"[MR] ERROR: Phenocode(s) not found in manifest: {', '.join(missing)}"
        if DEBUG: print(msg)
        return

    outdir = "local_dbs/FinnGen/"
    demo_outdir = "demo_outputs"
    os.makedirs(demo_outdir, exist_ok=True)

    # Ensure PLINK for clumping (no-op if already present)
    try:
        genal.install_plink()
    except Exception:
        pass

    # Fetch & load
    try:
        exp_file = _fetch_finngen_gwas(AVAIL_GWAS, exposure_gwas, outdir)
        out_file = _fetch_finngen_gwas(AVAIL_GWAS, outcome_gwas, outdir)
        exp_df_full = pd.read_csv(exp_file, sep="\t", compression="gzip", low_memory=False)
        out_df      = pd.read_csv(out_file,  sep="\t", compression="gzip", low_memory=False)
    except Exception as e:
        msg = f"[MR] ERROR fetching/reading GWAS: {e}"
        if DEBUG: print(msg)
        return

    # Optionally subset exposure to cis region around gene (±20 kb)
    exp_df = exp_df_full
    cis_mode = False
    cis_note = ""
    if gene_name:
        print('[MR] gene_name:', gene_name)
        genes_avail = list(gene_objs.keys())
        if gene_name not in genes_avail:
            print(f"[MR] Available genes for cis: {', '.join(genes_avail)}")
            raise ValueError(f"Gene '{gene_name}' not found in state['gene_entities']")
        else:
            gobj = gene_objs[gene_name]
            try:
                chrom, gstart, gend, strand = gobj.get_location()
            except Exception:
                chrom, gstart, gend, strand = None, None, None, None
            if (chrom is None) or (gstart is None) or (gend is None):
                if DEBUG: print(f"[MR] WARNING: gene '{gene_name}' missing hg38 coordinates; running genome-wide instruments.")
            else:
                try:
                    chrom_int = int(chrom)  # FinnGen uses integer #chrom on GRCh38
                    lo = max(0, int(gstart) - 20_000)
                    hi = int(gend) + 20_000
                    exp_df_cis = exp_df_full[(exp_df_full["#chrom"] == chrom_int) &
                                             (exp_df_full["pos"] >= lo) &
                                             (exp_df_full["pos"] <= hi)]
                    if exp_df_cis.empty:
                        if DEBUG:
                            print(f"[MR] WARNING: no exposure variants in ±20 kb for {gene_name} ({chrom_int}:{gstart}-{gend}); running genome-wide instruments.")
                    else:
                        print('Using cis variants:', exp_df_cis.shape[0])
                        exp_df = exp_df_cis
                        cis_mode = True
                        cis_note = f"As requested, limited the exposure instruments to variants within ±20 kb of {gene_name} (chr{chrom}:{lo}-{hi}). *Important* this changes interetation, where result reflects the effect of BMI changes driven by this gene region and not the genome-wide BMI effect, conclusion can be made that is only gene specific & careful as # variants might affect power."
                except Exception as e:
                    if DEBUG: print(f"[MR] WARNING: could not subset to cis for {gene_name}: {e}; running genome-wide instruments.")

    # Geno + preprocess (GRCh38 EUR)
    try:
        EXP = genal.Geno(exp_df, CHR="#chrom", POS="pos", EA="alt", NEA="ref",
                         BETA="beta", SE="sebeta", P="pval", EAF="af_alt", keep_columns=False)
        EXP.preprocess_data(preprocessing="Fill_delete", reference_panel="eur_38")
        EXP_clumped = EXP.clump(p1=5e-8, r2=0.1, kb=250, reference_panel="eur_38")
        if EXP_clumped is None:
            msg = f"[MR]{' [cis]' if cis_mode else ''} No instruments after clumping for '{exposure_gwas}'."
            message_state = f"Tried to run 2-sample MR with exposure '{exposure_gwas}' and outcome '{outcome_gwas}'" + (f" in cis mode around gene '{gene_name}'." if cis_mode else ".")+ "But no genome-wide significant instruments (p<5e-8) were found after clumping (r2=0.1, kb=250). Please try a different exposure."
            # add note to state
            state["text_notes"].append(message_state)
            if DEBUG: print(msg)
            return

        OUT = genal.Geno(out_df, CHR="#chrom", POS="pos", EA="alt", NEA="ref",
                         BETA="beta", SE="sebeta", P="pval", EAF="af_alt", keep_columns=False)
        OUT.preprocess_data(preprocessing="Fill_delete", reference_panel="eur_38")
    except Exception as e:
        msg = f"[MR]{' [cis]' if cis_mode else ''} ERROR during Geno/preprocess: {e}"
        if DEBUG: print(msg)
        return

    # Harmonize & MR
    try:
        EXP_clumped.query_outcome(OUT, proxy=False)
        res = EXP_clumped.MR(
            action=2,
            methods=["IVW", "IVW-FE", "WM", "Egger"],
            exposure_name=exposure_gwas,
            outcome_name=outcome_gwas,
            heterogeneity=True,
            odds=True,
        )
    except Exception as e:
        msg = f"[MR]{' [cis]' if cis_mode else ''} ERROR during MR: {e}"
        if DEBUG: print(msg)
        return

    # Descriptive single-paragraph summary (IDs, phenotypes, Ns, IVW+sig, heterogeneity, WM, Egger)
    try:
        summary = build_mr_summary_paragraph(AVAIL_GWAS, exposure_gwas, outcome_gwas, res, cohort_label="FinnGen R12", cis_mode=cis_mode, cis_comment=cis_note)
        state["text_notes"].append(summary)
    except Exception as e:
        msg = f"[MR]{' [cis]' if cis_mode else ''} ERROR summarizing: {e}"
        if DEBUG: print(msg)

    # Save structured payload for UI
    state.setdefault("two_sample_mr_ui", {})
    key = f"{exposure_gwas}_{outcome_gwas}" + (f"__cis_{gene_name}" if cis_mode else "")

    state["two_sample_mr_ui"][key] = {
        "cohort": "FinnGen R12",
        "summary_text": summary if 'summary' in locals() else None,
        "results_table": res,  # pandas DataFrame, as you requested
        "cis": cis_mode,
        "gene": gene_name if cis_mode else None,
    }




NODES: tuple[Node, ...] = (
    Node(
        name="twosample_mr_agent",
        entry_point=twosample_mr_agent,
        description=(
            "Run a 2-sample Mendelian Randomization (MR) analysis. "
            "Requires an exposure trait and an outcome trait (by FinnGen GWAS ID). "
            "Requires  gencode_gene_node to have been run first."
            "This analysis tests whether genetic variants *associated with the exposure "
            "causally influence the outcome* (very important which one is exposure and which one is outcome), helping to distinguish correlation from causation. "
            "Time consuming, can not run more that 2 combinations per session. "
            "Powerful for trait association, recommend running if user has any questions about available traits or connection needs to be done. "
            "TO RUN: istead of just `twosample_mr_agent` return as ('twosample_mr_agent-exposure_ID-outcome_ID') *with dashes*. "
            "*important*: if gene of interest is specified, make sure you ask for 'twosample_mr_agent-exposure_ID-outcome_ID-gene_name'."
            "Note that order is very important, first is *exposure*, second is outcome. "
            "Reverse order will change meaning of the test"
            "Only allowed one gene at a time, maximum 2 gene-commands per session. "
            "Note that IDs must be exact matched to the IDs shown below or it will fail. "
            "Available FinnGen phenotypes (ID: description): "
            + id_desc_str
        ),
        dependencies=('gencode_gene_node'),
    ),
)