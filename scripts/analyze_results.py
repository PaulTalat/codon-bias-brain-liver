import pandas as pd
from scipy import stats
import json

# 1. Load the Data
print("Loading results.csv...")
df = pd.read_csv('data/results.csv')

# 2. Separate the Groups
brain_data = df[df['group'] == 'Brain']['cai']
liver_data = df[df['group'] == 'Liver']['cai']

# 3. Calculate Descriptive Statistics
print("\n--- Descriptive Statistics ---")
print(f"Brain (n={len(brain_data)}): Mean CAI = {brain_data.mean():.3f}, Variance = {brain_data.var():.5f}")
print(f"Liver (n={len(liver_data)}): Mean CAI = {liver_data.mean():.3f}, Variance = {liver_data.var():.5f}")

# 4. Perform Statistical Testing (Mann-Whitney U Test)
# We use Mann-Whitney because sample sizes are small (n<30) and may not be perfectly normally distributed.
stat, p_value = stats.mannwhitneyu(brain_data, liver_data, alternative='two-sided')

print("\n--- Statistical Analysis ---")
print(f"Mann-Whitney U Statistic: {stat}")
print(f"P-value: {p_value:.5f}")

if p_value < 0.05:
    print("Result: SIGNIFICANT difference in codon bias between Brain and Liver genes.")
else:
    print("Result: NO SIGNIFICANT difference in codon bias between Brain and Liver genes.")
    
# 5. Extract top RSCU differences
# This block calculates the average RSCU for every codon in both groups to find the biggest differences.
print("\n--- Top RSCU Divergences ---")
brain_rscu_totals = {}
liver_rscu_totals = {}

for index, row in df.iterrows():
    rscu_dict = json.loads(row['rscu_json'])
    group = row['group']
    for codon, val in rscu_dict.items():
        if group == 'Brain':
            brain_rscu_totals.setdefault(codon, []).append(val)
        else:
            liver_rscu_totals.setdefault(codon, []).append(val)

differences = []
for codon in brain_rscu_totals.keys():
    brain_avg = sum(brain_rscu_totals[codon]) / len(brain_rscu_totals[codon])
    liver_avg = sum(liver_rscu_totals[codon]) / len(liver_rscu_totals[codon])
    diff = abs(brain_avg - liver_avg)
    differences.append((codon, diff, brain_avg, liver_avg))

differences.sort(key=lambda x: x[1], reverse=True)

for i in range(3):
    codon, diff, b_avg, l_avg = differences[i]
    print(f"{i+1}. {codon}: Diff={diff:.3f} (Brain Avg: {b_avg:.3f}, Liver Avg: {l_avg:.3f})")
    
# 6. Statistical Significance of Top RSCU Differences
print("\n--- RSCU Significance Testing ---")
# Testing our headline finding: Arginine codon divergence
top_codons_to_test = ['AGA', 'CGG', 'CTC']

for codon in top_codons_to_test:
    b_rscu = brain_rscu_totals[codon]
    l_rscu = liver_rscu_totals[codon]
    
    # Run the Mann-Whitney test specifically on this codon's arrays
    stat, p_val = stats.mannwhitneyu(b_rscu, l_rscu, alternative='two-sided')
    
    sig_status = "SIGNIFICANT" if p_val < 0.05 else "NOT significant"
    print(f"Codon {codon}: p-value = {p_val:.5f} ({sig_status})")