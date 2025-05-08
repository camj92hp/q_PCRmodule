## 📊 q_Pcr_Module


Thi is a Python-based module designed to streamline the analysis of quantitative PCR (qPCR) data using the ΔΔCt method. It allows researchers to clean and process raw qPCR outputs, normalize gene expression levels to reference genes and baseline samples, and visualize the resulting fold changes.

This tool is particularly useful for experiments comparing gene expression across treatments, timepoints, or genotypes, such as in immunology, cancer, or developmental biology research.


## 🧬 Features

- ✅ **Multi-file Support**: Load and combine multiple qPCR data files at once.
- 🧹 **Automated Data Cleaning**: Standardizes sample names, drops unnecessary columns, and removes invalid entries.
- 🧪 **ΔCt Calculation**: Computes ΔCt using reference genes like GAPDH.
- 🔁 **ΔΔCt Analysis**:
  - Treatment effects across timepoints
  - Gene mutation effects
  - Custom baseline normalization (e.g., control at 0 HR)
- 📊 **Visualization**: Bar plots showing log fold changes across genes, cell lines, and timepoints.
- 🔄 **Transposed Output**: Generates a wide-format version of the results for external statistical tools (e.g., Prism, R).
- 📁 **Excel Export**: Save results directly to `.xlsx` files for sharing and record keeping.
- 🧠 **Customizable Workflow**: Adaptable for different gene panels, timepoints, or naming conventions.


## ⚙️ Requirements
  This module is built for Python 3 and relies on the following libraries:

   🐍 Dependencies
  - pandas – DataFrame handling and file I/O
  - numpy – Array operation
  - seaborn – Data visualization
  - scipy – Statistical testing (e.g., t-tests)
  - matplotlib – Plot rendering and customization

  Make sure to install Python 3.7 or higher for full compatibility.

  You can install all required packages using pip:
    
        pip install -r requirements.txt

## 🚀 Usage
This toolkit is designed to make qPCR analysis as plug-and-play as possible.

### 🧪 Step 1: Prepare Your Data
Download the sample_data.xlsx file from this repository.

Replace the example content with your own qPCR output (exported from software like Bio-Rad or Thermo Fisher).

Make sure your file includes the columns:
    
    Sample, Target, Cq, and Biological Set Name.



### 📂 Step 2: Load and Clean the Data
````python
    
    from qpcr_analysis import load_file, data_cleaning

    # Load your file
    df = load_file("sample_data/sample_data.csv")

    # Define expected sample names
    samples = ["0 HR", "0.5 HR", "1 HR", "2 HRS", "4 HRS", "6 HRS"]

    # Clean and standardize
    clean_df = data_cleaning(samples, df)
````
### 🧬 Step 3: Run Analysis
Choose the analysis type that fits your experimental design.
All analysis functions let you specify your reference gene (ref) — default is "GAPDH".

#### ▶️ Treatment Effect
Compare timepoints or treatment conditions within the same cell line.

```python
    
    from qpcr_analysis import all_deltas_treatment_effect

    results = all_deltas_treatment_effect(clean_df, ref="GAPDH")  # Change ref gene if needed
 ```   
  - Input: One cell line over time or conditions
  - Reference: Timepoint 0 (via internal ΔCt)
  - Common Use: How does treatment change expression over time?

#### 🧬 Gene Mutation Effect
Compare the same timepoint across different genotypes or cell lines.

```python
    
    from qpcr_analysis import all_deltas_gene_effect
    results = all_deltas_gene_effect(clean_df, ref_cells='WT', ref="GAPDH")
```
  - Input: Multiple cell lines at same timepoint
  - Reference: A control cell line (e.g., 'WT')
  - Common Use: What is the effect of a mutation on expression at a fixed time?

#### 🕒 Baseline-Controlled Comparison
Normalize all expression values to a single baseline condition (e.g., WT at 0 HR), allowing comparison across both treatments and genotypes.

```python
  
    from qpcr_analysis import all_deltas_ctl_baseline

    results = all_deltas_ctl_baseline(clean_df, ref_cells='WT', ref_treatment='0 HR', ref="GAPDH")
```
  - Input: Multiple cell lines and timepoints
  - Reference: A fixed baseline (e.g., WT at 0 HR)
  - Common Use: How do different treatments and genetic backgrounds deviate from a single defined control

### 💾 Step 4: Save Results
```python
    from qpcr_analysis import saveToExcel
    saveToExcel(results, "results/qpcr_analysis_results.xlsx")
```
### 📊 Step 5: Visualize***
```python
    
    from qpcr_analysis import graph_qpcr
    graph_qpcr(results)
```


## 🧠 Core Functions

Below are the core functions included in this toolkit:

### 📂 File Handling & Cleaning

- `load_file(file_path)`
  - Loads a CSV file into a DataFrame. Assumes `index_col=0`.

- `data_cleaning(sample_name_list, *args)`
  - Concatenates and cleans multiple qPCR dataframes. Standardizes sample names, drops unused columns, and filters out NA values.

- `data_cleaning_single(df, drop, sample_name_list)`
  - Cleans a single dataframe with a custom list of columns to drop.

---

### 🧪 ΔCt & ΔΔCt Calculations

- `delta_ct(df, gene, cell, ref='GAPDH', sample='0 HR')`
  - Calculates ΔCt = GOI Cq − reference gene Cq for a given cell line and timepoint.

- `delta_delta_ct_treatment_effect(df, gene, cell, ref_cell, ref='GAPDH')`
  - Computes ΔΔCt to evaluate treatment effects across timepoints within a cell line.

- `delta_delta_ct_gene_effect(df, gene, cell, ref_cell='WT', treatment='6 HRS', ref='GAPDH')`
  - Compares gene expression between different cell types at the same treatment timepoint.

- `delta_delta_ct_ctl_baseline(df, gene, cell, ref_cell='WT', ref_treatment='0 HR', treatment='6 HRS', ref='GAPDH')`
  - Uses a fixed baseline (e.g., WT at 0 HR) to evaluate combined gene and treatment effects.

---

### 🔁 Batch ΔΔCt Wrappers

- `all_deltas_treatment_effect(df, ref='GAPDH')`
  - Runs treatment-effect ΔΔCt for all genes and cell types. Allows reference gene customization.

- `all_deltas_gene_effect(df, ref_cells='WT', ref='GAPDH')`
  - Compares all gene effects across cell lines and timepoints. Choose your baseline cell line and reference gene.

- `all_deltas_ctl_baseline(df, ref_cells='WT', ref_treatment='0 HR', ref='GAPDH')`
  - Combines gene and treatment effect using a fixed baseline control. Fully customizable.

---

### 📊 Visualization & Output

- `graph_qpcr(df)`
  - Creates bar plots of log fold changes for each gene across samples and cell types.

- `transposed_values(df)`
  - Reshapes data for external analysis (e.g., GraphPad Prism) with one row per sample group.

- `saveToExcel(df, filename)`
  - Exports any result DataFrame to an Excel file.


### 📈 Visualization
This toolkit includes built-in plotting functionality using Seaborn and Matplotlib.

🔬 graph_qpcr(df)
Generates a bar plot of log₂ fold change (log_change) values for each gene, grouped by sample and biological set (e.g., genotype or treatment group).

✅ Features:
-One subplot per gene
-Error bars showing standard deviation
-Grouped bars by biological condition
-Automatically adjusts layout for multiple genes

📍 Example:
```python

    from qpcr_analysis import graph_qpcr

    graph_qpcr(results)
```
This will display a figure like the one below:


### 🔄 Transposed Output Format
The transposed_values() function helps you quickly convert your results into a wide format that is ideal for external tools like GraphPad Prism or Excel-based statistics.

🧾 transposed_values(df)
Creates a row for each gene–sample–cell combination, and unpacks the log_change values into separate columns (replicates side-by-side), so you can easily copy and paste them into Prism.

📍 Example Output:
```python
      Target	Sample	Biological Set Name	log_change_1	log_change_2	log_change_3
      IFNβ	0 HR	WT	0.95	1.03	1.01
      IFNβ	6 HRS	KO	2.11	2.05	2.18
      TNFα	0 HR	WT	1.02	0.98	1.01
```
✅ This format is Prism-ready: simply copy each block of rows into a grouped bar graph or repeated measures layout.

🔧 How to Use:
```python

    from qpcr_analysis import transposed_values
    
    transposed_df = transposed_values(results)
```
You can then export the output:

```python
   
    from qpcr_analysis import saveToExcel
    
    saveToExcel(transposed_df, "results/prism_ready_output.xlsx")
```

### Saving Results

After analyzing and processing the data with the tool, users can easily save the results to an Excel file for further analysis or sharing. This functionality is built into the code through the `saveToExcel` function.

#### How to Save Your Results:

1. After performing the analysis, you can call the `saveToExcel` function.
2. Simply pass the DataFrame you wish to save along with the desired filename as arguments.
   
   ```python
   saveToExcel(df, "analysis_results.xlsx")


This will create an Excel file named analysis_results.xlsx containing the results in the current working directory.

If you want to save the results with a different file name, just provide a new name in the function:

```python
saveToExcel(df, "new_filename.xlsx")
```
The saved file will contain the processed results in a tabular format that can be easily imported into other tools like GraphPad Prism or Excel for further statistical analysis and visualization.
