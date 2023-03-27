import streamlit as st


def intro_page():
    st.markdown(
        """
        # Welcome to the RDKit

        RDKit is the industry standard, open-source package for cheminformatics.
        If you're interested in programmatically handling molecules, you are
        going to want to learn RDKit. This webpage itself uses RDKit for all of
        its chemical logic.

        You'll find that RDKit can be a little bit of a sprawling mess. One the
        one hand I hate that all the useful bits of RDKit are scattered across a
        massive code-base and difficult to search documentation. On the other, I
        love that it's powerful, fast, feature-rich, and easy to install and use
        once you feel comfortable with it.

        So it's time to jump in. I will introduce you to the most important
        parts of RDKit as well as let you in on some tips and tricks for using
        this powerful toolkit.

        Navigate to the next page using the side-bar.
        """
    )
