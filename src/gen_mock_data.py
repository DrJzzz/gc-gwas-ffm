import pandas as pd
import random

# Generate random phenotype dataset
num_samples = 100
phenotype_data = {
    "id_sample": [f"Sample{i+1}" for i in range(num_samples)],
    "phenotype": [round(random.uniform(0, 1), 2) for _ in range(num_samples)],
    "sex": [random.choice(["M", "F"]) for _ in range(num_samples)],
    "age": [random.randint(20, 60) for _ in range(num_samples)],
}

# Save to CSV
phenotypes = pd.DataFrame(phenotype_data)
phenotypes.to_csv("test_phenotypes.csv", index=False)

print("Phenotype dataset saved as 'test_phenotypes.csv'")
print(phenotypes)

# Use sample IDs from the phenotype file
sample_ids = phenotypes["id_sample"].tolist()

# Create the VCF file content
vcf_header = f"""##fileformat=VCFv4.2
##source=TestDataset
##contig=<ID=1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{'\t'.join(sample_ids)}
"""

vcf_body = ""
for variant_idx in range(1000):  # 100 variants
    pos = random.randint(1, 1000000)
    ref = random.choice(["A", "C", "G", "T"])
    alt = random.choice([a for a in ["A", "C", "G", "T"] if a != ref])
    genotypes = "\t".join(random.choice(["0/0", "0/1", "1/1"]) for _ in sample_ids)
    vcf_body += f"1\t{pos}\t.\t{ref}\t{alt}\t100\tPASS\t.\tGT\t{genotypes}\n"

# Write the VCF file
with open("test_genotypes.vcf", "w") as vcf_file:
    vcf_file.write(vcf_header + vcf_body)

print("Synthetic VCF file 'test_genotypes.vcf' created!")