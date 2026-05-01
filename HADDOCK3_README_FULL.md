# HADDOCK3 Installation Guide

## Overview

This guide describes a **reliable and reproducible method** to install **HADDOCK3 (v2026.3.0)** 

Two critical issues must be addressed:

- **Environment contamination (`PYTHONPATH`)**
- **Avoiding source compilation of dependencies**

This procedure has been tested and works on:

- CentOS 7
- Ubuntu 22.04.5 LTS
- AlmaLinux 9.6  
- Rocky Linux 9.7

---

## Requirements & Common Pitfalls

### 1. Environment Contamination (`PYTHONPATH`)

Some HPC environments (e.g., Amber installations) modify the `PYTHONPATH` variable. This can cause Python to load **incompatible packages** from external installations.

This leads to:
- installation failures
- runtime errors
- inconsistent environments

#### 🔍 Check your environment:

```bash
echo $PYTHONPATH
```

If you see paths like:

```
/usr/local/bio/amber20/...
```

then your environment is contaminated.

### Clean it:

```bash
unset PYTHONPATH
```

>  **Important:** A polluted `PYTHONPATH` was a major source of installation failures.

---

## 2. Installation steps

### Step 2.1 — Create a clean conda environment

Use Python ≥ 3.10 (recommended: 3.11)

```bash
conda create -n haddock3 python=3.11 pandas numpy scipy -c conda-forge
```

---

### Step 2.2 — Activate the environment

```bash
conda activate haddock3
```

---

### Step 2.3 — Ensure clean environment


```bash
unset PYTHONPATH
```

---

### Step 2.4 — Install HADDOCK3

```bash
pip install haddock3
```

---

## Step 3 — Verify installation

```bash
haddock3 -h
```

or:

```bash
python -c "import haddock"
```

---
---
---


## If the current version of HADDOCK3 is incompatible with the ligand search, follow the instructions below
### Install HADDOCK3 enviroment 2026.03.0

To ensure reproducibility, we provide:

---

### Recreate the environment

```bash
conda env create -f haddock3_environment_for_sharing.yml
conda activate haddock3
```

---

### Best Practices

- Always use a **clean conda environment**
- Do **not** mix system Python with conda environments
- Avoid globally setting `PYTHONPATH`
- Keep Amber (or similar tools) isolated using modules or separate scripts
- Prefer conda for scientific dependencies and pip only for HADDOCK3

---

