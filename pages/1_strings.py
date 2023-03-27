import streamlit as st

from pages.strings.intro import intro_page
from pages.strings.smarts import smarts_1_page, smarts_2_page
from pages.strings.smiles import smiles_page
from pages.strings.smirks import smirks_1_page

subpages = [
    ("0 - intro", intro_page),
    ("1 - SMILES", smiles_page),
    ("2 - SMARTS Part 1", smarts_1_page),
    ("3 - SMARTS Part 2", smarts_2_page),
    ("4 - SMIRKS Part 1", smirks_1_page),
]
st.sidebar.markdown("## sub-topics")
subpage = st.sidebar.selectbox(
    "select page", subpages, format_func=lambda x: x[0]
)
st.sidebar.markdown("---")
subpage[1]()
