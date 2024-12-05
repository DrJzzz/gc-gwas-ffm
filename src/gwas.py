### Forma esperada del dataset:
# 1. Informacion del fenotipo: <id_sample>, <phenotype>, <edad>, <sexo>
# 2. Informacion de los genotipos: se espera que sea un archivo VCF

import pandas as pd
import allel
import statsmodels.api as sm
import numpy as np

# import phenotype data:
phenotypes = pd.read_csv('test_phenotypes.csv')


# asegurarse de que esten normalizados los datos
if phenotypes['phenotype'].min() < 0 or phenotypes['phenotype'].max() > 1:
    raise ValueError("El rango de los valores fenotipicos no esta entre 0 y 1.")


# import genotype data:

g_dataset = allel.read_vcf('test_genotypes.vcf')

sample_ids = g_dataset['samples']
positions = g_dataset['variants/POS']
chromosomes = g_dataset['variants/CHROM']
genotypes = g_dataset['calldata/GT']


## Asegurarse de que los datos esten alineados
phenotype_ids = phenotypes['id_sample'].values
if not all(phenotype_ids == sample_ids):
    raise ValueError("Los ids de los samples no coinciden entre el archivo de genotipos y el de fenotipos.")

phenotypes_filtered = phenotypes[phenotypes['id_sample'].isin(sample_ids)]

aligned_indices = [i for i, sample in enumerate(sample_ids) if sample in phenotypes_filtered['id_sample'].values]
genotypes_aligned = genotypes.take(aligned_indices, axis=0)

###

phenotype_values = phenotypes['phenotype'].values

genotype_data = genotypes_aligned.to_n_alt()

## Tambien se examinan los valores de edad y sexo
covariates = phenotypes_filtered[['edad', 'sexo']].values