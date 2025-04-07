# **📊 SFPPySt: Interactive Tools for Food Contact Assessment using Streamlit**

<div align="center">

| ![Generative Simulation](https://github.com/ovitrac/SFPPylite/raw/main/extra/assets/logo.png) | This Streamlit app extends the [Generative Simulation](https://github.com/ovitrac/generativeSimulation) initiative | Say it.<br />Simulate it with AI. |
| ------------------------------------------------------------ | ------------------------------------------------------------ | --------------------------------- |

</div>

---

> 🍏⏩🍎 ***SFPPy**: A Python Framework for Food Contact Compliance & Risk Assessment*  

### 🚀 SFPPy-St: Streamlit-powered SFPPy Tools

**Fast, interactive graphical interfaces for SFPPy** — built with [**Streamlit**](https://streamlit.io), deployed server-side or locally:

- 🧪 **Substance query** and **toxicological summary** via `migrantToxtree`
- 🧬 Display of EU, US, and CN regulatory data for ~1300 known substances
- 🧾 Live PubChem queries (structure, synonyms, InChIKey, SMILES, etc.)
- 💡 Designed to complement [**SFPPyLite**](https://github.com/ovitrac/SFPPylite) with full toxicological support (e.g. ToxTree)

[![SFPPy_st App](https://img.shields.io/badge/SFPPy_st-Search%20Substance-4CAF50?logo=streamlit&logoColor=white)](https://sfppyst.streamlit.app/) (*click the badge to launch the app.)*

---

### 🎯 Key Differences from SFPPyLite

| Feature             | **SFPPy-St (Streamlit)**                                     | **SFPPyLite (JupyterLite)**                                  |
| ------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| **PubChem Support** | Server-side Python with full `loadpubchem` support for not documented substances | In-browser Python (Pyodide/WASM), `loadpubchem` cannot cache new substances |
| **ToxTree Support** | ✅ Full support for Cramer classification, alerts, TTC if installed locally. <br />⚠️ Support for cached substances only with the on line server (no JVM on Streamlit servers). | ⚠️ Support for cached substances (Java dependency)            |
| **Interactivity**   | GUI forms, dropdowns, live search                            | Jupyter notebook UI                                          |
| **Performance**     | Native Python speed                                          | ~2x slower (WebAssembly overhead)                            |
| **Persistence**     | ✅ if user-local, ❌ with Streamlit                            | In-browser (IndexedDB) with 📥📤                               |
| **Deployable**      | Yes (locally or on a web server)                             | Yes (static website)                                         |
| **Mobile ready**    | ✅ (Streamlit is responsive)                                  | ✅ (Lite runs in-browser)                                     |
| **Launching time**  | <5 s if user-local, < 2 min with Streamlit                   | <10 s                                                        |

---

### 🛠 Features in SFPPy-St

- 🔍 **Search substances** using name, CAS, or synonyms
- 📊 Visualize **Cramer class, TTC values, and alerts**
- 🧾 Access **regulatory frameworks** from EU, FDA, and China
- 🖼 View **PubChem structure images** with size control and auto-rotation
- 📈 Ready for deployment with `streamlit run app.py`

---

### 🧪 Launch SFPPy-St (Local)

1. **Install Streamlit** (if not already):
```bash
conda install -c conda-forge streamlit
```

2. **Run the app**:
```bash
streamlit run search_substance.py
```

---

### 🍏⏩🍎 Related Projects

[![SFPPy](https://img.shields.io/badge/SFPPy-%F0%9F%8D%8F%E2%8F%A9%F0%9F%8D%8E_PARENT%20PROJECT-4CAF50?style=for-the-badge&logo=python)](https://github.com/ovitrac/SFPPy)
[![SFPPyLite](https://img.shields.io/badge/SFPPyLite-Launch%20in%20Browser-blueviolet?logo=jupyter&style=for-the-badge)](https://github.com/ovitrac/SFPPyLite)

- **SFPPy**: Full-featured desktop framework with simulation, fitting, compliance rules
- **SFPPyLite**: Lightweight browser version of SFPPy for notebooks, without JVM dependencies
- **SFPPy-St**: This repo — GUI-friendly version of SFPPy for toxicological exploration

---

### 💡 Why Streamlit?

- ⚡ Rapid GUI prototyping
- 📦 No frontend code needed
- 🌍 Ideal for demonstration, education, or deployment
- 📜 Fully compatible with SFPPy codebase

---

### 📬 Feedback

💬 Found a bug or want to suggest a tool? [Open an issue](https://github.com/ovitrac/SFPPy/issues) or email the author.

---

### 📁 Project

Core tools and databases are located in `patankar/`. Do not rename this folder.

```
SFPPySt/
├── search_substance.py      # Main app file (Streamlit)
├── patankar/                # Optional: grouped apps
    ├── private/             # private modules and databases shipped with SFPPy
    ├── cache.PubChem/       # Pubchem cache 
    ├── cache.ToxTree/       # ToxTree cache
    └── loadpubchem.py       # substance manager connected to PubChem
    └── geometry.py          # packaging geometry module
    └── food.py              # food module
    └── layer.py             # polymer/material module
    └── migration.py         # mass transfer solver
    └── property.py          # transfport/thermodynamic module
    └── useroverride.py      # user override
├── requirements.txt         # (for pip)
├── environment.yml          # (for conda)
├── README.md                # this file
└── LICENSE.md               # MIT License
```

---

### 🧰 Powered by

- [Streamlit](https://streamlit.io)
- [SFPPy](https://github.com/ovitrac/SFPPy)
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
- [Pillow](https://pillow.readthedocs.io/) for structure image manipulation