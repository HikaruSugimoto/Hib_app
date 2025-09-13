
import os, io, re, glob, base64
from typing import Dict, List, Optional
from PIL import Image
import streamlit as st
from streamlit.components.v1 import html as st_html
from pathlib import Path
import zipfile
import sys
from typing import Tuple, Dict, Optional
import pandas as pd
import numpy as np

st.set_page_config(page_title="A multi-tissue, multi-omics atlas of hibernation", layout="wide")
if(os.path.isfile('demo.zip')):
    os.remove('demo.zip')
with zipfile.ZipFile('demo.zip', 'x') as csv_zip:
    csv_zip.writestr("demo.csv",
                    pd.read_csv("demo.csv").to_csv(index=False))    
with open("demo.zip", "rb") as file:
    #st.sidebar.download_button(label = "Download demo data",data = file,file_name = "demo.zip")
    zip_data = file.read()
    b64 = base64.b64encode(zip_data).decode()
    zip_filename = 'demo.zip'
    href = f'<a href="data:application/zip;base64,{b64}" download="{zip_filename}">Download demo data</a>'
    st.sidebar.markdown(href, unsafe_allow_html=True)
if(os.path.isfile('demo.zip')):
    os.remove('demo.zip')

st.title("A multi-tissue, multi-omics atlas of hibernation")
st.markdown(
    """
This app compares your uploaded omics results with a curated multi-omics database of hibernation- or hypometabolism-associated changes 
(15,771 alterations across transcriptome/proteome/miRNA/metabolome from 285 studies). This outputs overlapping molecules and categorizes them by directional agreement.
    """
)

st.write('This app accepts data in the following format:')
image = Image.open('input data.png')
st.image(image, caption='',use_container_width=True)

#prior knowledge
curated_raw = None
load_errs = []
for enc in ["cp932", "shift_jis", "utf-8", "utf-8-sig"]:
    try:
        curated_raw = pd.read_csv("previous.csv", encoding=enc)
        break
    except Exception as e:
        load_errs.append(f"{enc}: {e}")
if "Organ" in curated_raw.columns:
    curated_raw["Organ"] = curated_raw["Organ"].replace({
        "Muscle (SOL)": "Muscle","Muscle (EDL)": "Muscle",
        "Muscle (Plantaris)": "Muscle","Muscle (GAS)": "Muscle",
        "Myotube": "Muscle",'Jejunum': 'Intestine','Cecal epithelium': 'Intestine',
        'Hypothalamus': 'Brain','Pituitary': 'Brain','Cortex': 'Brain','Midbrain': 'Brain',
        'Forebrain': 'Brain','Medulla': 'Brain','Hippocampus': 'Brain','Platelet': 'Blood','Red blood cell':'Blood',
        'Adipocyte': 'WAT'
    })
organ = st.sidebar.selectbox("Organ", ['Brain','Retina', 'Heart','Lung', 
   'Liver', 'Kidney','Adrenal gland','Spleen', 'Intestine', "Ovary",
   'Wing', 'Bone', 'WAT','BAT',
   'Muscle','Aorta','Blood'], index=0)
assay_type = st.sidebar.selectbox("Type", ["mRNA", "protein", "miRNA", "metabolite"], index=0)
curated = curated_raw.copy()
if "Type" in curated.columns:
    curated = curated[curated["Type"].astype(str).str.lower() == assay_type.lower()]
if "Organ" in curated.columns:
    curated = curated[curated["Organ"].astype(str).str.lower() == organ.lower()]

if curated.empty:
    st.warning("Curated database is empty after filtering. Try different filters.")
if "Molecule" in curated.columns:
    curated["Molecule"] = curated["Molecule"].astype(str).str.strip()
if "Prior_Up_Down" in curated.columns:
    curated["Prior_Up_Down"] = curated["Prior_Up_Down"].astype(str).str.strip().str.lower().map({"up":"Up", "down":"Down"})
else:
    st.error("'UP_Down' column is missing in previous.csv; cannot compute overlaps.")
    st.stop()
curated_up = curated[curated["Prior_Up_Down"] == "Up"].copy()
curated_down = curated[curated["Prior_Up_Down"] == "Down"].copy()


#upload data        
user_upload = st.sidebar.file_uploader("Upload your dataset (CSV)", type=["csv"], key="user_csv")
if user_upload is None:
    st.info("➡️ Upload your CSV in the sidebar to proceed.")
    st.stop()

user_df = None
user_errs = []
for enc in ["utf-8", "utf-8-sig", "cp932", "shift_jis", "shift_jisx0213"]:
    try:
        user_df = pd.read_csv(user_upload, encoding=enc)
        break
    except Exception as e:
        user_errs.append(f"{enc}: {e}")
        user_upload.seek(0)

if user_df is None or user_df.empty:
    st.error("Could not read your CSV. Tried encodings: " + ", ".join([e.split(":")[0] for e in user_errs]))
    st.stop()

st.subheader("Uploaded dataset")
st.dataframe(user_df.head(20), use_container_width=True)

user_std = pd.DataFrame()
user_std["Molecule"] = user_df["Molecule"].astype(str).str.strip()
raw = user_df["Up_Down"].astype(str).str.strip().str.lower()
ups = {"up", "upregulated", "+", "increase", "increased"}
downs = {"down", "downregulated", "-", "decrease", "decreased"}
user_std["User_UP_Down"] = pd.Series(pd.NA, index=user_df.index, dtype="string")
user_std.loc[raw.isin(ups), "User_UP_Down"] = "Up"
user_std.loc[raw.isin(downs), "User_UP_Down"] = "Down"

user_std = user_std.dropna(subset=["User_UP_Down", "Molecule"]).drop_duplicates(subset=["Molecule"]).copy()
if user_std.empty:
    st.warning("No rows classified as Up or Down. Check your column selection and thresholds.")
    st.stop()

# -----------------------
up_up = pd.merge(user_std[user_std["User_UP_Down"] == "Up"], curated_up, on="Molecule", how="inner", suffixes=("_user", "_prior"))
up_down = pd.merge(user_std[user_std["User_UP_Down"] == "Up"], curated_down, on="Molecule", how="inner", suffixes=("_user", "_prior"))
down_down = pd.merge(user_std[user_std["User_UP_Down"] == "Down"], curated_down, on="Molecule", how="inner", suffixes=("_user", "_prior"))
down_up = pd.merge(user_std[user_std["User_UP_Down"] == "Down"], curated_up, on="Molecule", how="inner", suffixes=("_user", "_prior"))

concordant = pd.concat([up_up, down_down], ignore_index=True)
discordant = pd.concat([up_down, down_up], ignore_index=True)

prior_up_set = set(curated_up["Molecule"].dropna().unique())
prior_down_set = set(curated_down["Molecule"].dropna().unique())
both_prior = prior_up_set & prior_down_set
user_all_set = set(user_std["Molecule"].unique())
controversial = pd.DataFrame({"Molecule": sorted(user_all_set & both_prior)})
any_overlap = pd.DataFrame({"Molecule": sorted(user_all_set & (prior_up_set | prior_down_set))})

# -----------------------
# Summary metrics
# -----------------------
n_all = user_std["Molecule"].nunique()
concordant_set = set(concordant["Molecule"]) if not concordant.empty else set()
discordant_set = set(discordant["Molecule"]) if not discordant.empty else set()
controversial_set = set(controversial["Molecule"]) if not controversial.empty else set()

colA, colB, colC, colD = st.columns(4)
with colA:
    st.metric("All unique molecules", n_all)
with colB:
    st.metric("Concordant (same direction)", len(concordant_set))
with colC:
    st.metric("Discordant (opposite direction)", len(discordant_set))
with colD:
    st.metric("Controversial in prior", len(controversial_set))

st.markdown("---")

# -----------------------
# Tabs and downloads
# -----------------------
tabs = st.tabs([
    "Concordant (Same)",
    "Discordant (Opposite)",
    "Controversial (Prior Up & Down)",
    "Any Overlap",
    "Raw: Up-Up", "Raw: Up-Down", "Raw: Down-Down", "Raw: Down-Up"
])

(tab_conc, tab_disc, tab_contr, tab_any, tab_uupu, tab_updo, tab_dndn, tab_dnup) = tabs

with tab_conc:
    st.dataframe(concordant, use_container_width=True)
    #if not concordant.empty:
    #    st.download_button("Download concordant.csv", concordant.to_csv(index=False).encode("utf-8"), file_name="concordant.csv", mime="text/csv")

with tab_disc:
    st.dataframe(discordant, use_container_width=True)
    #if not discordant.empty:
    #    st.download_button("Download discordant.csv", discordant.to_csv(index=False).encode("utf-8"), file_name="discordant.csv", mime="text/csv")

with tab_contr:
    st.dataframe(controversial, use_container_width=True)
    #if not controversial.empty:
    #    st.download_button("Download controversial.csv", controversial.to_csv(index=False).encode("utf-8"), file_name="controversial.csv", mime="text/csv")

with tab_any:
    st.dataframe(any_overlap, use_container_width=True)
    #if not any_overlap.empty:
    #    st.download_button("Download any_overlap.csv", any_overlap.to_csv(index=False).encode("utf-8"), file_name="any_overlap.csv", mime="text/csv")

with tab_uupu:
    st.dataframe(up_up, use_container_width=True)
    #if not up_up.empty:
    #    st.download_button("Download up_up.csv", up_up.to_csv(index=False).encode("utf-8"), file_name="up_up.csv", mime="text/csv")

with tab_updo:
    st.dataframe(up_down, use_container_width=True)
    #if not up_down.empty:
    #    st.download_button("Download up_down.csv", up_down.to_csv(index=False).encode("utf-8"), file_name="up_down.csv", mime="text/csv")

with tab_dndn:
    st.dataframe(down_down, use_container_width=True)
    #if not down_down.empty:
    #    st.download_button("Download down_down.csv", down_down.to_csv(index=False).encode("utf-8"), file_name="down_down.csv", mime="text/csv")

with tab_dnup:
    st.dataframe(down_up, use_container_width=True)
    #if not down_up.empty:
    #    st.download_button("Download down_up.csv", down_up.to_csv(index=False).encode("utf-8"), file_name="down_up.csv", mime="text/csv")

st.markdown("---")

# Optional: offer a ZIP with all results
zbuf = io.BytesIO()
with zipfile.ZipFile(zbuf, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
    if not concordant.empty: zf.writestr("concordant.csv", concordant.to_csv(index=False))
    if not discordant.empty: zf.writestr("discordant.csv", discordant.to_csv(index=False))
    if not controversial.empty: zf.writestr("controversial.csv", controversial.to_csv(index=False))
    if not any_overlap.empty: zf.writestr("any_overlap.csv", any_overlap.to_csv(index=False))
    if not up_up.empty: zf.writestr("up_up.csv", up_up.to_csv(index=False))
    if not up_down.empty: zf.writestr("up_down.csv", up_down.to_csv(index=False))
    if not down_down.empty: zf.writestr("down_down.csv", down_down.to_csv(index=False))
    if not down_up.empty: zf.writestr("down_up.csv", down_up.to_csv(index=False))
    if not curated.empty: zf.writestr("previous.csv", curated.to_csv(index=False))
zbuf.seek(0)
st.download_button(
    label="📦 Download all results (results.zip)",
    data=zbuf.getvalue(),
    file_name="results.zip",
    mime="application/zip",
)

