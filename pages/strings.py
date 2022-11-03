import streamlit as st
from PIL import Image

from pages.strings_pages import (
    intro,
    smiles,
)

subpages = [
    ('0 - intro', intro.page),
    ('1 - SMILES', smiles.page)
]

st.sidebar.markdown('## sub-topics')
subpage = st.sidebar.selectbox(
    'select page',
    subpages,
    format_func=lambda x: x[0]
)
st.sidebar.markdown("---")
subpage[1]()
