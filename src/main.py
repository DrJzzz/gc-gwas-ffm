import pandas as pd
import allel
import statsmodels.api as sm
import numpy as np
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

# Cargar datasets 
phenotype_file = "test_phenotypes.csv"
phenotypes = pd.read_csv(phenotype_file)

vcf_file = "test_genotypes.vcf"
callset = allel.read_vcf(vcf_file)

# extraer info del vcf
sample_ids = callset['samples']
positions = callset['variants/POS']
chromosomes = callset['variants/CHROM']
genotypes = allel.GenotypeArray(callset['calldata/GT'])

# checar alineacion
aligned_indices = [i for i, sample in enumerate(sample_ids) if sample in phenotypes['id_sample'].values]
genotypes_aligned = genotypes.take(aligned_indices, axis=0)
phenotypes_aligned = phenotypes[phenotypes['id_sample'].isin(sample_ids)]
phenotype_values = phenotypes_aligned['phenotype'].values

# Print dataset summary
print(f"Number of samples: {len(aligned_indices)}")
print(f"Number of variants: {len(positions)}")


### GWAS

# Flatten genotype data (convert to allele counts: 0, 1, or 2)
genotype_data = genotypes_aligned.to_n_alt()

# Initialize lists to store GWAS results
gwas_results = []

# Run linear regression for each SNP
for snp_idx in range(genotype_data.shape[1]):
    # Genotype for current SNP
    snp = genotype_data[:, snp_idx]
    
    # Prepare design matrix (include intercept)
    X = sm.add_constant(snp)
    y = phenotype_values
    
    # Perform linear regression
    model = sm.OLS(y, X).fit()
    
    # Collect results
    gwas_results.append({
        'Chromosome': chromosomes[snp_idx],
        'Position': positions[snp_idx],
        'Effect_Size': model.params[1],
        'P_Value': model.pvalues[1]
    })

# Convert results to DataFrame
gwas_df = pd.DataFrame(gwas_results)

# Display first few results
print(gwas_df.head())


### Bonferroni correction

# Apply FDR correction
gwas_df['Adjusted_P_Value'] = multipletests(gwas_df['P_Value'], method='fdr_bh')[1]

# Filter significant SNPs
significant_snps = gwas_df[gwas_df['Adjusted_P_Value'] < 0.05]
print("Significant SNPs:")
print(significant_snps)

### PLOTS

# Add -log10(P_Value) for visualization
gwas_df['-log10(P_Value)'] = -np.log10(gwas_df['P_Value'])

# Plot Manhattan plot
plt.figure(figsize=(12, 6))
plt.scatter(gwas_df['Position'], gwas_df['-log10(P_Value)'], c='blue', s=2)
plt.axhline(y=-np.log10(5e-8), color='red', linestyle='--')  # Genome-wide significance threshold
plt.xlabel("Genomic Position")
plt.ylabel("-log10(P-value)")
plt.title("Manhattan Plot")
plt.show()

# QQ plot
observed = -np.log10(sorted(gwas_df['P_Value']))
expected = -np.log10(np.linspace(1 / len(gwas_df), 1, len(gwas_df)))

plt.figure(figsize=(6, 6))
plt.scatter(expected, observed, c='blue', s=2)
plt.plot([0, max(expected)], [0, max(expected)], color='red', linestyle='--')
plt.xlabel("Expected -log10(P-value)")
plt.ylabel("Observed -log10(P-value)")
plt.title("QQ Plot")
plt.show()


### SAVE
gwas_df.to_csv("gwas_results.csv", index=False)
print("GWAS results saved to 'gwas_results.csv'")