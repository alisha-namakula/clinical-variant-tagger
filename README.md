
# clinical-variant-tagger

An R pipeline to tag pathogenic *TP53* missense mutations using their biophysical amino acid properties.

# 🧬 *TP53* Clinical Variant Pathogenicity Tagger Pipeline Using Biophysical Features

## 📌 Summary

I built an R pipeline to analyse and tag protein variants from clinical data of patient missense mutations in the TP53 gene. I used this to suggest clinical conclusions for patients. I wanted to improve my R skills over the holidays and fuel a deep interest in bioinformatics.

---

## ⚙️ How It Works

### 1. Data Cleaning

First, I took the data and removed all columns I would not need. I kept the codon, as well as the name and single-letter code of the original wild-type (`WT_AA`) and mutant amino acids (`Mutant_AA`).

### 2. Left Join

I completed the dataset by attaching the patient missense mutations to their corresponding amino acid properties. This included molecular weight (mass), hydrophobicity (oiliness), and pI (charge).

### 3. Feature Engineering

I took the complete clean dataset and mutated it to calculate the percentage change between `Mutant_AA` and `WT_AA`. This allowed me to quantify the shift in oiliness and molecular mass to determine their clinical significance.

The "charge" columns for `WT_AA` and `Mutant_AA` were not transformed into percentage changes. Charge swaps of any degree are almost always detrimental, so I did not create a rule to account for a margin of error.

### 4. Clinical Tagger Rules

I used the `case_when()` function to impose conditional rules to classify the data into distinct categories to aid in determining clinical significance.

These were the rules used to classify them accordingly:

* `is.na(Mol_Wt_Pct_change) | is.na(Hydro_ph7_Pct_change)` = **"CHECK FOR FRAMESHIFT"**
* `WT_AA Charge != Mut_AA Charge` = **PATHOGENIC** (`!=` means "not equal to")
* `Hydrophobicity % Change < -30` = **PATHOGENIC**
* `Hydrophobicity % Change > 30` = **PATHOGENIC**
* `Molecular Weight % Change < -30` = **PATHOGENIC**
* `Molecular Weight % Change > 30` = **PATHOGENIC**
* `Remaining` = **BENIGN**

### 5. Data Visualisation

I used `ggplot2` to generate 3 types of graphs: bar charts, histograms, and a scatter plot to display pathogenicity trends and biophysical changes across variants.

---

## 📊 Graph Details

### **Graph 1: Bar Chart**

* **Title:** *"Distribution of Pathogenicity Across Patient Missense Mutations"*

### **Graph 2a: Histogram (Hydrophobicity)**

* **Title:** *"Normal Distribution of Hydrophobicity Changes Across Missense Mutations"*

### **Graph 2b: Histogram (Molecular Weight)**

* **Title:** *"Most Prevalent Molecular Weight Change Across Missense Mutations"*

### **Graph 3: Scatter Plot (Molecular Weight vs. Hydrophobicity)**

* **Title:** *"Molecular Weight and Hydrophobicity Percentage Changes in Missense Mutations"*
* **Subtitle:** *"Sub-Contributing Factors to Pathogenicity"*

---

## 📈 Visual Data Trends

* **Bar Chart:** Shows how many patient mutations ended up tagged as Pathogenic vs. Benign vs. Check for Frameshift. Many missense mutations were tagged as benign.
* **Hydrophobicity Histogram:** Shows that most mutations have a near 0% change in oiliness. The extreme cases were more spread out.
* **Molecular Weight Histogram:** Shows that most amino acid swaps were close to their original size or heaviness as well, having a near 0% change.
* **Scatter Plot:** Shows how molecular weight and hydrophobicity changes look together for each individual mutation. It further visualises how these two factor changes huddle close to 0 in most missense mutations.

---

## 🧬 Biological & Clinical Insights

Most of the swaps in the missense mutations were benign.

Pathogenicity in these *TP53* mutations may not be mostly attributed to changes in hydrophobicity or molecular weight, as the majority of amino acid swaps kept these two factors quite similar — meaning they may not cause severe damage to the overall folding of the *TP53* protein. Instead, these variants may mainly be tagged as pathogenic because charge acts as the primary driver of disease in *TP53* missense mutations.

---

## 📁 References

* **Patient Mutation Dataset:** UMD *TP53* Mutation Database (Clinical cohort of human *TP53* missense variants).
* **Amino Acid Biophysical Properties:** Sigma-Aldrich / IUPAC Amino Acid Reference Data (Molecular Weight, Hydropathy Index at pH 7, and Isoelectric Point pI).

---

## 🛠️ Tools Used

* **Language:** R
* **Packages:** `dplyr`, `ggplot2`

---
