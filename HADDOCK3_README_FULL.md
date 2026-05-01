# 🧬 HADDOCK3 Installation Guide (CentOS 7)

## 📌 Overview

This guide describes a **reliable and reproducible method** to install **HADDOCK3 (v2026.3.0)** on **CentOS 7** and similar HPC environments.

Two critical issues must be addressed:

- ⚠️ **Environment contamination (`PYTHONPATH`)**
- ⚠️ **Avoiding source compilation of dependencies**

This procedure has been tested and works on:

- CentOS 7  
- AlmaLinux 9.6  
- Rocky Linux 9.7  

---

## ⚠️ Requirements & Common Pitfalls

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

#### 🧹 Clean it:

```bash
unset PYTHONPATH
```

> 💡 **Important:** A polluted `PYTHONPATH` was a major source of installation failures.

---

### 2. Avoid Source Compilation (CentOS 7 Limitation)

Running:

```bash
pip install haddock3
```

may trigger **source compilation** of packages like:

- pandas  
- numpy  
- scipy  

This often fails on CentOS 7 due to:

- outdated compilers  
- missing headers (e.g., `stdatomic.h`)  
- incompatible system toolchains  

#### ✅ Solution

Avoid compilation entirely by installing scientific dependencies via **conda** (precompiled binaries).

---

## 🚀 Installation Steps (Recommended)

### Step 1 — Create a clean conda environment

Use Python ≥ 3.10 (recommended: 3.11)

```bash
conda create -n haddock3 python=3.11 pandas numpy scipy -c conda-forge
```

---

### Step 2 — Activate the environment

```bash
conda activate haddock3
```

---

### Step 3 — Ensure clean environment

```bash
unset PYTHONPATH
```

---

### Step 4 — Install HADDOCK3

```bash
pip install haddock3
```

✔️ This works because:

- dependencies are already installed via conda  
- pip does **not** need to compile anything  
- installation is fast and stable  

---

### Step 5 — Verify installation

```bash
haddock3 -h
```

or:

```bash
python -c "import haddock"
```

---

## 🔁 Reproducing the Environment

To ensure reproducibility, we provide:

```
haddock3_environment_for_sharing.yml
```

This file was generated using:

```bash
conda env export --no-builds > haddock3_environment_for_sharing.yml
```

---

### Why `--no-builds`?

- ✅ Removes system-specific build strings  
- ✅ Improves portability across clusters  
- ✅ Keeps package versions consistent  

---

### 📦 Recreate the environment

```bash
conda env create -f haddock3_environment_for_sharing.yml
conda activate haddock3
```

---

## 🧠 Best Practices

- Always use a **clean conda environment**
- Do **not** mix system Python with conda environments
- Avoid globally setting `PYTHONPATH`
- Keep Amber (or similar tools) isolated using modules or separate scripts
- Prefer conda for scientific dependencies and pip only for HADDOCK3

---

## ✔️ Summary

A successful HADDOCK3 installation on CentOS 7 requires:

- Python ≥ 3.10 (recommended: 3.11)  
- Clean environment (`PYTHONPATH` unset)  
- Conda-installed scientific stack (no compilation)  
- pip installation of HADDOCK3  

---

## 💬 Notes

This guide reflects a **tested and working configuration** on HPC systems where:

- system compilers are outdated  
- mixed software environments are common  
- reproducibility is important  

---

