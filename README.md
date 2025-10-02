# Single-Cell Lung Transplant Project

## 📌 Overview

This project analyzes **ischemia-reperfusion responses in human lung
transplants** at the **single-cell resolution** using **scRNA-seq**. By
profiling donor and recipient lung tissues collected before and after
reperfusion, we aim to uncover cell-type-specific inflammatory responses
that drive lung injury during transplantation.

## 🧬 Objectives

-   Process and analyze single-cell RNA sequencing (scRNA-seq) data from
    donor and recipient lung samples.\
-   Identify **differential gene expression** across timepoints (pre-
    vs. post-reperfusion) and between donor vs. recipient cells.\
-   Characterize key immune cell populations (e.g., **macrophages,
    monocytes, CD8+ T cells**) driving inflammatory responses.\
-   Map cytokine/chemokine signaling pathways (e.g., **IL-1β, IL-6,
    CXCL8**) involved in lung injury.\
-   Explore **cell-cell interaction networks** to understand
    donor--recipient immune crosstalk.

## 🛠️ Methods

-   **Preprocessing & QC**: Standard scRNA-seq filtering, normalization,
    clustering, and annotation.\
-   **Differential Expression**: Compare pre- vs. post-reperfusion, and
    donor vs. recipient cells.\
-   **Pathway Analysis**: Identify enriched inflammatory signaling
    pathways.\
-   **Cell-Cell Interaction**: Predict ligand--receptor interactions
    between donor and recipient immune cells.

## 📊 Expected Outcomes

-   A cell atlas of donor and recipient immune responses during
    ischemia-reperfusion.\
-   Insights into **molecular drivers of graft injury** at single-cell
    resolution.\
-   Potential therapeutic targets for reducing lung injury in
    transplantation.

## 💻 Usage

1.  Clone the repository

    ``` bash
    git clone https://github.com/yourusername/single-cell-lung-transplant.git
    cd single-cell-lung-transplant
    ```

2.  Install dependencies (e.g., Seurat for R or Scanpy for Python).

    ``` bash
    pip install scanpy anndata matplotlib seaborn
    ```

3.  Run the analysis scripts provided in the `scripts/` directory.

## 📂 Data

Raw and processed data are available upon request or from
controlled-access repositories (e.g., GEO, dbGaP).

## 📜 License

This project is licensed under the MIT License. See `LICENSE` for
details.

## ✨ Acknowledgements

Based on research by Wong et al., *American Journal of Transplantation
(2024)* and associated collaborators at the University of Toronto.
