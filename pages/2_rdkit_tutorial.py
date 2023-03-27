import streamlit as st

from pages.rdkit_tutorial.introduction import intro_page
from pages.rdkit_tutorial.mol import mol_page

subpages = [
    ("0 - intro", intro_page),
    ("1 - mol", mol_page),
]
st.sidebar.markdown("## sub-topics")
subpage = st.sidebar.selectbox(
    "select page", subpages, format_func=lambda x: x[0]
)
st.sidebar.markdown("---")
if subpage is not None:
    subpage[1]()
