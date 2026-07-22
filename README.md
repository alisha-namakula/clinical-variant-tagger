# clinical-variant-tagger
An R pipeline to tag pathogenic TP53 missense mutations using biophysical amino acid properties.

# 🧬 Clinical Variant Pathogenicity Tagger & Biophysical Pipeline

## 📌 Summary 
I built an RStudio pipeline to tag clinical data of patient missense mutations from the gene **TP53**. I wanted to improve my R skills over the holidays and fuel a deep interest in bioinformatics.

---

## ⚙️ How It Works 
### 1. Data Cleaning
First, I took the data and removed all columns I would not need. I kept the codon, as well as the name and single-letter code of the original (WT_AA) and mutant amino acids (Mutant_AA).

### 2. Left Join
I completed the dataset by attaching the patient missense mutations to their corresponding amino acid properties. This included molecular weight (mass), hydrophobicity (oiliness), and pI (charge).

### 3. Feature Engineering
I took the complete clean dataset and mutated it to calculate the percentage change between `Mutant_AA` and `WT_AA`. This allowed me to quantify the shift between oiliness and molecular mass to determine their clinical significance. 

The "charge" `WT_AA` and `Mutant_AA` columns were not transformed into percentage changes. Charge swaps of any degree are almost always detrimental, so I did not create a rule to account for a margin of error.

### 4. Logic Tagger
I used the `case_when()` function to impose conditional rules to classify the data into 4 distinct categories to aid in determining clinical significance.

These were the rules to classify them accordingly:
* `is.na(Mol_Wt_Pct_change) | is.na(Hydro_ph7_Pct_change)` = **"CHECK FOR FRAMESHIFT"**
* `WT_AA Charge != Mut_AA Charge` = **PATHOGENIC**
* `Hydrophobicity % Change < -30` = **PATHOGENIC**
* `Hydrophobicity % Change > 30` = **PATHOGENIC**
* `Molecular Weight % Change < -30` = **PATHOGENIC**
* `Molecular Weight % Change > 30` = **PATHOGENIC**
* `Remaining` = **BENIGN**

---

## 📊 Data Visualization

I used `ggplot2` to generate 3 types of graphs: bar charts, histograms, and a scatter plot.

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

## 💡 Main Takeaways from the Graphs
* **Bar Chart:** Shows how many patient mutations ended up tagged as Pathogenic vs. Benign vs. Frameshift.
* **Hydrophobicity Histogram:** Shows that most mutations have near 0% change in oiliness, but the extreme ones spread out.
* **Molecular Weight Histogram:** Shows that most amino acid swaps were also close to the same size or heaviness (0% change).
* **Scatter Plot:** Shows how molecular weight and hydrophobicity changes look together for each mutation.

---

## 📁 References
* **Patient Mutation Dataset:** UMD TP53 Mutation Database (Clinical cohort of human *TP53* missense variants).
* **Amino Acid Biophysical Properties:** Sigma-Aldrich / IUPAC Amino Acid Reference Data (Molecular Weight, Hydropathy Index at pH 7, and Isoelectric Point pI).

---

## 🛠️ Tools Used
* **Language:** R
* **Packages:** `dplyr`, `ggplot2`
