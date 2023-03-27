import streamlit as st

from streamlit_utils import (
    st_check_smirks,
    st_display_reaction,
    st_run_reaction,
)


def smirks_1_page():
    st.markdown(
        """
        # SMIRKS

        As far as I can tell, unlike SMILES or SMARTS, SMIRKS doesn't actually
        stand for anything. That means I get to make something up! How about
        SMarts-Imitating ReaKtion String. So far we've covered string
        representations of molecules (SMILES) and substructures (SMARTS) but my
        favorite cheminformatic string representation is definitely the SMIRKS.

        ## Atom mapping

        Let's start by looking at one of the simplest reaction in organic
        chemistry, the reaction of methyl bromide with sodium ethoxide.
        """
    )

    reactants = ["CBr", "[O-]CC"]
    products = ["COCC", "[Br-]"]
    st_display_reaction(reactants=reactants, products=products)

    st.markdown(
        """
        In this classic SN2 reaction, a `C-Br` bond is being broken and a new
        `O-C` bond is being formed. As it turns out, SMIRKS has a way to encode
        this very reaction in a string format!

        SMIRKS is very similar to SMARTS. Actually, the only difference between
        SMIRKS and SMARTS are:
        * you must have exactly 2 `>` characters in your string
            * the format is `[reactants]>[reagents]]>[products]
        * you must give the atom mapping (reactant atoms to product atoms)
            * the format is `:#` where `#` is the atom number
        """
    )

    smk = "[C:1][Br:2].[O-:3][C:4][C:5]>>[C:1][O-0:3][C:4][C:5].[Br-:2]"
    reactants = ["CBr", "[O-]CC"]
    st_run_reaction(smk, reactants)

    line_too_long = (
        "as we see in the atom numbers, `[C:1]` comes from the methyl bromide w"
        "hile the rest of the atoms come from ethoxide"
    )
    st.markdown(
        f"""
        Let's break down this SMIRKS:
        * `[C:1][Br:2]` represents methyl bromide
            * the carbon is now atom 1 and the bromine is atom 2
        * `[O-:3][C:4][C:5]` represents the ethoxide
            * the atoms in ethoxide are now assigned an atom number as well
        * `[C:1][O:3][C:4][C:5]` represents the product, ethyl methyl ether
            * {line_too_long}
        * `[Br-:2]` represents bromide, the by-product of this reaction

        What makes SMIRKS powerful is that it is built off of SMARTS. It's meant
        to be generalizable, and again, capable of expressing our "chemical
        intuition" of what should and shouldn't react.

        For example, we actually don't need to include all the atoms on the
        ethoxide. We don't need both atoms `4` and `5`, we only actually need
        `4`.
        """
    )

    smk = "[C:1][Br:2].[O-:3][C:4]>>[C:1][O-0:3][C:4].[Br-:2]"
    st_run_reaction(smk, reactants)

    st.markdown(
        """
        This is because atom `5` is not involved in the reaction. We should also
        rewrite this SMIRKS to ignore the by-product bromide ion.
        """
    )

    smk = "[C:1]Br.[O-:3][C:4]>>[C:1][O-0:3][C:4]"
    st_run_reaction(smk, reactants)

    st.markdown(
        """
        And there we have it, a SMIRKS string that does a pretty good job at
        encoding a basic SN2 reaction.

        ### Your turn!

        I know that writing SMIRKS is somewhat more tedious than SMILES or
        SMARTS. But I would strongly encourage you to give these a serious try
        before reaching for the answer button.
        """
    )

    prompt = "Reaction: Etherification"
    reactantss = [["CCl", "Oc1ccccc1"]]
    productss = [["COc1ccccc1"]]
    answer = "[C:1]Cl.[O:2][c:3]>>[C:1][O:2][c:3]"
    st_check_smirks(reactantss, productss, answer, prompt)

    prompt = "Reaction: Reductive amination"
    reactantss = [["CC=O", "NC1CC1"]]
    productss = [["CCNC1CC1"]]
    answer = "[C:1]=O.[N:2][C:3]>>[C:1][N:2][C:3]"
    st_check_smirks(reactantss, productss, answer, prompt, idx_offset=1)

    st.markdown(
        """
        ## Making our SMIRKS more smart

        What makes SMIRKS really powerful is the power of SMARTS that it
        inherits. For example, what if we wanted a SMIRKS that could capture
        reaction of alkoxides with methyl chloride and bromide?
        """
    )

    smk = "[C:1][Cl,Br].[O:2][C:3]>>[C:1][O:2][C:3]"
    reactants = ["CCl", "OC"]
    st_run_reaction(smk, reactants)
    st.markdown("---")
    reactants = ["CBr", "OC"]
    st_run_reaction(smk, reactants, display_smk=False)

    st.markdown(
        """
        As we see, the key component here is the `[Cl,Br]`. This gives the
        SMIRKS some flexibility to select either methyl chloride or bromide.

        What if we wanted to go the other way and add some restrictions to our
        SMIRKS? For example this SMIRKS allows tert-butyl alcohol to react but
        we know that SN2 is less likely for hindered nucleophiles.
        """
    )

    smk = "[C:1][Cl,Br].[O:2][C:3]>>[C:1][O:2][C:3]"
    reactants = ["CCBr", "OC"]
    st_run_reaction(smk, reactants)
    st.markdown("---")
    reactants = ["CCBr", "OC(C)(C)C"]
    st_run_reaction(smk, reactants, display_smk=False)

    st.markdown(
        """
        What if we wanted to exclude t-butyl alcohols? We could add in some
        SMARTS properties to `[C:3]`.
        """
    )

    smk = "[C:1][Cl,Br].[O:2][C;!H0:3]>>[C:1][O:2][C:3]"
    reactants = ["CCBr", "OC"]
    st_run_reaction(smk, reactants)
    st.markdown("---")
    reactants = ["CCBr", "OC(C)(C)C"]
    st_run_reaction(smk, reactants, display_smk=False)

    st.markdown(
        """
        Specifying the carbon attached to the alcohol as `[C;!H0:3]` excludes
        the t-butyl alcohol as that particular carbon has no hydrogens. However,
        it's not that simple, what about phenols? Our SMIRKS excludes phenols
        but we know that phenols should be able to react.
        """
    )

    reactants = ["CCBr", "Oc1ccccc1"]
    st_run_reaction(smk, reactants)

    line_too_long_1 = (
        "We are specifying the carbon attached to the oxygen as `C`. This will "
        "exclude aromatic carbons, like the carbons in the phenol ring."
    )
    line_too_long_2 = (
        "We are requiring the carbon attached to the oxygen to have 1, 2, or 3 "
        "hydrogens but the carbon attached to the phenol oxygen has 0."
    )
    st.markdown(
        f"""
        There are two problems.
        * {line_too_long_1}
        * {line_too_long_2}

        This can be fixed by using `#6` in place of `C` and using recursion to
        aliphatic carbon logic from aromatic carbon logic.
        """
    )

    smk = "[C:1][Cl,Br].[O:2][#6;$([C!H0]),$(c):3]>>[C:1][O:2][#6:3]"
    reactants = ["CCBr", "OC"]
    st_run_reaction(smk, reactants)
    st.markdown("---")
    reactants = ["CCBr", "OC(C)(C)C"]
    st_run_reaction(smk, reactants, display_smk=False)
    st.markdown("---")
    reactants = ["CCBr", "Oc1ccccc1"]
    st_run_reaction(smk, reactants, display_smk=False)

    st.markdown(
        """
        Hopefully you see now why SMIRKS is my personal favorite. Using the
        flexibility and specificity of SMARTS, we're able to encode the logic of
        a chemical reaction into a sequence of letters, numbers, and symbols.

        ### Your turn!
        """
    )

    prompt = "Reaction: SNAr of aryl fluorides."
    reactantss = [
        ["c1cccnc1F", "N1CCCC1"],
        ["c1ccccc1F", "N1CCCC1"],
    ]
    productss = [["c1cccnc1N1CCCC1"], [""]]
    answer = "F[c:1][n:2].[N:3]>>[N:3][c:1][n:2]"
    st_check_smirks(reactantss, productss, answer, prompt, idx_offset=2)

    prompt = "Reaction: Boc deprotection."
    reactantss = [
        ["C1CCN1C(=O)OC(C)(C)C"],
    ]
    productss = [
        ["C1CCN1"],
    ]
    answer = "[#7:1]C(=O)OC(C)(C)C>>[#7:1]"
    st_check_smirks(reactantss, productss, answer, prompt, idx_offset=3)

    prompt = "Reaction: Reduction of ketones."
    reactantss = [
        ["C1C(=O)CCC1"],
        ["C1C(=O)NCC1"],
        ["C1C(=O)OCC1"],
    ]
    productss = [
        ["C1C(O)CCC1"],
        [""],
        [""],
    ]
    answer = "[C;$(C(C)C):1]=[O:2]>>[C:1][O:2]"
    st_check_smirks(reactantss, productss, answer, prompt, idx_offset=4)

    st.markdown(
        """
        And there we have it, a primer on SMIRKS. Once you have enough practice
        with SMIRKS, many doors will open. SMIRKS is used for many, many tasks
        in cheminformatics:
        * exchanging between different tautomers
        * virtually deprotecting molecules
        * building virtual compound libraries
        """
    )
