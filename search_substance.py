import streamlit as st
import pandas as pd
from PIL import Image
from patankar.loadpubchem import migrantToxtree, migrant



def auto_rotate_if_tall(image):
    """
    Rotate the image by 90 degrees if it is taller than it is wide.
    """
    if image.height > image.width:
        return image.rotate(90, expand=True)
    return image


def show_attr(label, value):
    """generic display"""
    if value is None or (isinstance(value, str) and value.strip() == ""):
        return

    # Handle CID with PubChem link
    if label == "CID":
        url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{value}"
        st.markdown(f"**{label}:** [{value}]({url})")
        return

    # Handle short list → inline
    if isinstance(value, list):
        if len(value) == 0:
            return
        elif len(value) < 3:
            st.markdown(f"**{label}:** {', '.join(str(v) for v in value)}")
            return
        else:
            df = pd.DataFrame({label: value})
            st.markdown(f"**{label}:**")
            st.dataframe(df, height=min(35 * len(df), 350))
            return

    # Handle dict → vertical key-value table
    if isinstance(value, dict):
        if not value:
            return
        df = pd.DataFrame(list(value.items()), columns=["Key", "Value"])
        st.markdown(f"**{label}:**")
        st.dataframe(df, height=min(35 * len(df), 350))
        return

    # Scalar fallback
    st.markdown(f"**{label}:** {m.dispmax(value)}")



def display_migrant(m):
    """Display a migrant or migrantToxtree instance in Streamlit."""


    def section(title): st.markdown(f"### {title}")

    # General Info
    section("🧬 General Information")
    show_attr("Compound", m.compound)
    show_attr("Synonyms", m.name)  # if m.name is a list

    # Show image
    if hasattr(m, "image") and m.image:
        img = auto_rotate_if_tall(m.image)
        st.image(img, caption=f"Structure of {m.compound}", width=450, output_format="PNG")


    # Handle name: dropdown for synonyms
    # if isinstance(m.name, list):
    #     st.markdown("**Synonyms (select one):**")
    #     selected_name = st.selectbox("Synonyms", options=m.name)
    #     st.markdown(f"**Selected name:** {selected_name}")
    # else:
    #     show_attr("Name", m.name)

    show_attr("CID", m.cid)
    show_attr("CAS", m.CAS)
    st.markdown("---")
    show_attr("M", f"{m.M} g/mol")
    show_attr("Formula", m.formula)
    show_attr("SMILES", getattr(m, "smiles", None))
    show_attr("InChiKey", getattr(m, "InChiKey", None))
    st.markdown("---")
    show_attr("logP", m.logP)
    show_attr("P' (calc)", m.polarityindex)

    # 🇪🇺 EC 10/2011
    if m.hasSML:
        st.markdown("<hr style='border: 1px solid #ccc; margin-top: 1em; margin-bottom: 1em;'>",unsafe_allow_html=True)
        section("🇪🇺 EC 10/2011")
        if m.hasannex1:
            show_attr("SML", f"{m.SML} [{m.SMLunit}]")
            if m.annex1["SMLTGroupFCMsubstances"]:
                show_attr("Group", f"{len(m.annex1['SMLTGroupFCMsubstances'])} substances")
            show_attr("Annex Name", m.annex1["name"])
            show_attr("Annex CAS", m.annex1["CAS"])
            show_attr("EC|FCM|Ref", f"{m.annex1['EC']}|{m.annex1['FCM']}|{m.annex1['Ref']}")
        else:
            show_attr("SML", f"{m.SML} [{getattr(m, 'SMLunit', '')}]")

    # 🇺🇸 FCN
    if m.hasFCN:
        st.markdown("<hr style='border: 1px solid #ccc; margin-top: 1em; margin-bottom: 1em;'>",unsafe_allow_html=True)
        section("🇺🇸 US FDA FCN")
        show_attr("FCM No", m.FCNNo)
        show_attr("Notifier", m.FCNnotifier)
        show_attr("Manufacturer", m.FCNmanufacturer)
        show_attr("Notification Date", m.FCNnotificationDate)
        if m.FCNmixture:
            show_attr("Mixture", f"{m.FCNnsubstances} substances")

    # 🇨🇳 GB9685-2016
    if m.hasFCA:
        st.markdown("<hr style='border: 1px solid #ccc; margin-top: 1em; margin-bottom: 1em;'>",unsafe_allow_html=True)
        section("🇨🇳 CN GB9685-2016")
        show_attr("FCA No", m.FCANo)
        show_attr("Authorized in", m.FCAgroups)
        if "plastics" in m.FCAgroups:
            show_attr("Polymers", m.FCApolymers)
            if m.FCACP0max: show_attr("🇨🇳 CP0 max", f"{m.FCACP0max} mg/kg")
            if m.FCASML: show_attr("🇨🇳 SML", f"{m.FCASML} mg/kg")
            if m.FCAQM: show_attr("🇨🇳 QM", f"{m.FCAQM} mg/kg")
            if m.FCADL: show_attr("🇨🇳 DL", f"{m.FCADL} mg/kg")

    # 🧪 ToxTree
    if isinstance(m, migrantToxtree) and m.compound not in (None, "", []):
        st.markdown("<hr style='border: 1px solid #ccc; margin-top: 1em; margin-bottom: 1em;'>",unsafe_allow_html=True)
        section("🧪 ToxTree Results")
        show_attr("IUPAC Traditional Name", m.ToxTree.get("IUPACTraditionalName"))
        show_attr("IUPAC Name", m.ToxTree.get("IUPACName"))
        show_attr("Cramer Class", m.CramerClass)
        show_attr("TTC", f"{m.TTC} {m.TTCunits}")
        show_attr("CF TTC", f"{m.CFTTC} {m.CFTTCunits}")

        # Show alerts only if non-empty
        alerts = m.showalerts if isinstance(m.showalerts, dict) else {}
        if any(alerts.values()):
            st.markdown("#### Alerts")
            st.json(alerts)

# %% main()
# Optional: set custom Streamlit page config
st.set_page_config(page_title="SFPPy Search Substance", layout="centered")

# --- Title and input ---
st.title("🧪 SFPPy — Search Substance")
st.markdown("Enter a **compound name**, **synonym**, or **CAS number** to fetch data from PubChem and evaluate toxicology.")

compound_name = st.text_input("Substance name or CAS", value="BPA")

if st.button("Search & Run ToxTree") and compound_name.strip():
    try:
        with st.spinner("Running ToxTree..."):
            m = migrant(compound_name, verbose=False)

        st.success(f"✅ Substance loaded: {m.compound}")
        mtox = m.promote()
        if not isinstance(mtox, migrantToxtree):
            st.warning("⚠️ No cached ToxTree data available. Toxicology section skipped.")
        st.divider()
        display_migrant(mtox)

    except Exception as e:
        st.error(f"❌ Error: {e}")
