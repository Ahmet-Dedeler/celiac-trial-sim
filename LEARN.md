# Learning Clinical Trial Simulation: Biostatistics & Power Analysis

This guide explains the biostatistical concepts behind the Celiac Trial Simulator.

---

## 1. The Core Scientific Problem

In celiac disease trials with gluten food challenges (e.g., TAK-101, latiglutenase, PRV-015), primary endpoints measure mucosal damage via duodenal biopsy:
- **Villous Height to Crypt Depth Ratio (Vh:Cd)**
- **Intraepithelial Lymphocyte (IEL) counts**

Several trials failed to meet primary statistical endpoints and concluded the drugs were ineffective. However, this simulation demonstrates that:
1. **Endpoint noise is not constant:** Standard deviation scales with mucosal injury (\(\text{SD} = 0.40 + 0.30 \times \text{injury}\)).
2. **Underpowered sample sizes:** With only 13 patients per arm, a drug would need an unrealistic 95% protection efficacy to achieve \(p < 0.05\) at 80% power. TAK-101 prevented 71% of injury yet appeared as a statistical failure due to small cohort sizes.

---

## 2. Simulation Methodology

- **Synthetic Patient Cohorts:** Generates Monte Carlo cohorts sampling baseline histology, gluten challenge intake, and biological response heterogeneity.
- **Heteroskedastic Noise:** Models non-linear variance expansion during inflammatory mucosal blunting.
- **Power Sweeps:** Performs thousands of simulated trials across varying cohort sizes (\(N=10\) to \(N=150\)) to construct statistical power curves.
