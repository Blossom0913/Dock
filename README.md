# Dock

A molecular docking study investigating potential ligand candidates against the 5O3L receptor protein. 

## 👥 Contributors

This project is a collaboration between:
- **Xutian Chen**
- **Yuanyuan Wei**

*UCLA*

## 🎯 Project Overview

This repository contains molecular docking analysis using the 5O3L receptor protein against multiple compound libraries to identify potential binding candidates. 

## 🧬 Receptor Protein

- **Target:** 5O3L. pdbqt
- **Format:** PDBQT (Protein Data Bank, Partial Charge (Q), & Atom Type (T))

## 📚 Ligand Libraries

The study screens ligand candidates from the following databases: 

| Database | Number of Molecules | Processing Time | Throughput |
|----------|-------------------|-----------------|------------|
| **Genesis** | 102,726 | 11 hours (RTX 2080 Ti) | 2.6 mol/sec |
| **MIPE5.0** | 2,490 | - | - |
| **Sytravon** | 44,953 | - | - |
| **NPACT** | 5,100 | - | - |
| **NPC-2019** | 2,678 | - | - |
| **Enamine** | 32,000,000 | - | - |

**Total Compounds Screened:** ~32.2 million molecules

## 🔬 Methodology

This project employs computational molecular docking techniques to:
1. Screen large chemical libraries against the 5O3L receptor
2. Identify potential binding candidates
3. Rank compounds based on predicted binding affinity
4. Facilitate drug discovery and lead optimization

## 💻 Hardware

- **GPU:** NVIDIA RTX 2080 Ti
- **Performance:** 2.6 molecules per second (Genesis library benchmark)

## 📁 Repository Structure

```
Dock/
├── receptor/
│   └── 5O3L.pdbqt          # Target receptor protein
├── ligands/                 # Ligand libraries
│   ├── genesis/
│   ├── mipe5.0/
│   ├── sytravon/
│   ├── npact/
│   ├── npc-2019/
│   └── enamine/
├── results/                 # Docking results and analysis
├── scripts/                 # Docking and analysis scripts
└── README.md
```

## 🚀 Getting Started

### Prerequisites

- AutoDock Vina or similar docking software
- Python 3.x
- CUDA-compatible GPU (recommended for large-scale screening)

### Installation

```bash
git clone https://github.com/Blossom0913/Dock.git
cd Dock
```

### Usage

*(Add specific instructions for running your docking simulations)*

```bash
# Example command
# vina --receptor receptor/5O3L.pdbqt --ligand ligands/example.pdbqt --out results/output.pdbqt
```

## 📊 Results

*(Add summary of key findings, top-ranked compounds, or visualization of results)*

## 🔗 References

### Databases
- **Genesis:** [Add link if available]
- **MIPE5.0:** Molecular Integrated Platform for Exploration version 5.0
- **Sytravon:** [Add link if available]
- **NPACT:** Naturally Occurring Plant-based Anti-cancer Compound-Activity-Target database
- **NPC-2019:** NCBI Natural Products Collection
- **Enamine:** [Enamine Chemical Database](https://enamine.net/)

### 5O3L Protein
- [PDB Entry for 5O3L](https://www.rcsb.org/structure/5O3L)

## 📝 License

*(Add your preferred license)*

## 📧 Contact

For questions or collaboration inquiries, please open an issue or contact the contributors. 

---

**Note:** This project is part of research conducted at UCLA. Please cite appropriately if using this work. 
