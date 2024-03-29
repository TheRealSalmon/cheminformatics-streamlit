import streamlit as st

from streamlit_utils import st_check_smarts, st_display_mol_with_smi


def smarts_1_page():
    st.markdown(
        """
        # SMARTS Part 1

        SMARTS are a SMiles ARbitrary Target Specification. SMARTS are an
        extension of SMILES that allows even more precise specification of
        chemical structure. While SMILES focuses on representing molecules,
        SMARTS focuses on representing substructures.

        In Part 1, we will introduce the basic syntax which is somewhat more
        complicated than SMILES. In Part 2, however, we will really get to see
        the power and flexibility of SMARTS with the introduction of recursive
        definitions.

        ## Why SMARTS? Substructure specification is hard.

        In general, SMILES is insufficient to specify a molecular substructure.
        In some cases, it is possible to specify the substructure you have in
        mind with just a SMILES string. For example, quaternary ammoniums can be
        entirely specified with just SMILES.
        """
    )

    smi = "C[N+](C)(C)C"
    smiles = ["CN", "CNC", "CN(C)C", "C[NH+](C)C"]
    st_display_mol_with_smi(smiles, substruct=smi, use_smiles=True)
    smiles = [
        "C[N+](C)(C)C",
        "C[N+](C)(C)C(C)(C)C",
        "C[N+](C)(CC1)C1",
        "C[N+](C)(C)c1ccccc1",
    ]
    st_display_mol_with_smi(
        smiles, substruct=smi, use_smiles=True, display_substruct=False
    )

    st.markdown(
        """
        As it turns out, the SMILES string `C[N+](C)(C)C` is entirely capable of
        differentiating quaternary ammoniums from primary, secondary, and
        tertiary amines.

        However, there are cases where SMILES is insufficient. Consider below
        the case of separating carboxylic acids from esters.
        """
    )

    smi = "O=C[OH]"
    smiles = ["O=CO", "CC(=O)O", "CC(=O)OC"]
    st_display_mol_with_smi(smiles, substruct=smi, use_smiles=True)

    st.markdown(
        """
        As we see, despite putting an `H` on the oxygen of the carboxylic acid,
        we're unable to separate acids from esters.
        """
    )

    st.info(
        """
        Despite being explicitly written into the SMILES as `[OH]`, RDKit does
        not explicitly match the `H` when checking substructures.
        """
    )

    st.markdown(
        """
        We *could* solve this with SMILES by using programming logic to require
        the substructure to match `O=C[OH]` *and* not match `O=COC`. This would
        filter out esters. However, how would we deal with carbamates? What
        about carbonates? Wouldn't we have to write a unique filtering logic for
        each substructure?
        """
    )
    smiles = ["CNC(=O)OC", "COC(=O)OC"]
    st_display_mol_with_smi(smiles, substruct=smi, use_smiles=True)

    st.markdown(
        """
        This is where SMARTS come to the rescue. SMARTS are similar to SMILES
        but they have the characteristics of being more flexible *and* more
        specific.
        """
    )

    smt = "O=C[OH]"
    smiles = ["O=CO", "CC(=O)O", "CC(=O)OC"]
    st_display_mol_with_smi(smiles, substruct=smt)
    smiles = ["CNC(=O)OC", "COC(=O)OC"]
    st_display_mol_with_smi(smiles, substruct=smt, display_substruct=False)

    st.markdown(
        """
        As we see here, SMARTS is more strict when it comes to counting the
        number of Hs on the oxygen, allowing us to separate carboxylic acids
        from esters.

        In the rest of this tutorial, we will introduce the SMARTS syntax,
        focusing on atom properties.

        ## Atoms

        In SMILES we were somewhat restricted in the atomic properties that
        could be specified. However, in SMARTS, the list is much longer:
        * already in SMILES:
            * charge
            * isotope
            * &#8203;# of hydrogens
            * chirality
            * aromatic/aliphatic
        * wildcard: `[A]`, `[a]`, or `*`
        * degree or # of connections: `[CX#]`
        * if atom is part of a ring: `[CR]` or `[C!R]`
        * the size of the ring the atom is in: `[Cr#]`

        ## Wildcards

        SMARTS allows us to include wildcards in our substructures. There are a
        few types of wildcards in SMARTS:
        * `A` will match aliphatic atoms.
        * `a` will match aromatic atoms.
        * `*` is a general wildcard, it will match *any* atom.
        """
    )

    smt = "A"
    smiles = ["CNC", "c1ocnc1", "C1NCCC1n1cccc1"]
    st_display_mol_with_smi(smiles, substruct=smt)

    smt = "a"
    st_display_mol_with_smi(smiles, substruct=smt)

    smt = "*"
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown(
        """
        SMARTS strictly separates aliphatic and aromatic atoms. For example,
        with carbon, `C` only matches aliphatic carbons while `c` matches
        aromatic carbons.
        """
    )

    smt = "C"
    st_display_mol_with_smi(smiles, substruct=smt)

    smt = "c"
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown(
        """
        However, what if we wanted a SMARTS atom that matched both aliphatic
        and aromatic carbons? In that case, we can specify the carbon by using
        its atomic number, `[#6]`.
        """
    )

    smt = "[#6]"
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown(
        """
        Lastly, SMARTS allows us to create our own wildcards by using SMARTS'
        most powerful feature, logic. For example, we may want a "middle-ground"
        of flexibility in-between specifying atomic number and using a wildcard.
        This is where the `,` character comes in.

        What if we want to match benzenes and pyridines but not anything else?
        """
    )

    smt = "[c,n]1ccccc1"
    smiles = ["c1ccccc1", "n1ccccc1", "[o+]1ccccc1", "p1ccccc1"]
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown(
        """
        Using the `,` keyword, we could tell the SMARTS to allow `c` *or* `n`
        only. However, this is just the beginning of SMARTS logic!

        ### Your turn!
        """
    )

    prompt = (
        "Separate the compounds containing nitrogens from the compounds that do"
        " not contain nitrogen."
    )
    smiles_in = ["CNCC", "c1ncccc1", "C1CCC1C(=O)NNC"]
    smiles_out = ["COCC", "c1occc1", "C1CCC1C(=O)OO"]
    answer = "[#7]"
    st_check_smarts(smiles_in, smiles_out, answer, prompt)

    prompt = (
        "Separate the 5-membered aromatic rings from the 6-membered aromatic ri"
        "ngs."
    )
    smiles_in = ["c1occc1", "c1[nH]ccc1", "c1sccc1"]
    smiles_out = ["c1ccccc1", "c1ncccc1", "c1[o+]cccc1"]
    answer = "a1aaaa1"
    st_check_smarts(smiles_in, smiles_out, answer, prompt, idx_offset=1)

    prompt = (
        "Separate the reactive epoxides and aziridines from the more stable cyc"
        "lopropanes."
    )
    smiles_in = ["c1occc1C1OC1", "C1CCC(N2)C12", "CC1N(CC(=O)O)C1C"]
    smiles_out = ["c1occc1C1CC1", "C1CCC(C2)C12", "CC1C(CC(=O)O)C1C"]
    answer = "C1[O,N]C1"
    st_check_smarts(smiles_in, smiles_out, answer, prompt, idx_offset=2)

    st.markdown(
        """
        ## Atomic degree

        Let's go back to the example of the quaternary ammonium above. In SMILES
        we can use `C[N+](C)(C)C`. However, this explicitly writes out each
        carbon. With SMARTS, we can specify the degree and the # of hydrogens on
        the `N` to write a more concise substructure string.
        """
    )

    smt = "[NX4H0+]"
    smiles = ["CN", "CNC", "CN(C)C", "C[NH+](C)C"]
    st_display_mol_with_smi(smiles, substruct=smt)
    smiles = [
        "C[N+](C)(C)C",
        "C[N+](C)(C)C(C)(C)C",
        "C[N+](C)(CC1)C1",
        "C[N+](C)(C)c1ccccc1",
    ]
    st_display_mol_with_smi(smiles, substruct=smt, display_substruct=False)

    st.markdown(
        """
        Breaking it down, the `X4` indicates that this `N` has `4` neighbors.
        The `H0` says there are no hydrogens on this `N` and the `+` gives us
        the charge.

        In another case, `X#` can be used to separate imine-like nitrogens from
        amine-like nitrogens.
        """
    )

    st.info(
        """
        The `X` keyword counts the number of hydrogen and non-hydrogen
        connections.
        """
    )

    smt = "[#7X2]"
    smiles = ["CNC", "c1ccccc1N", "CC=N", "c1nc[nH]c1"]
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown(
        """
        ### Your turn!

        Try to solve these problems using only the `X` keyword.
        """
    )

    prompt = "Separate the carbonyls from the ethers/alcohols."
    smiles_in = ["CC=O", "C1CNC(=O)C1", "c1c(=O)[nH]ccc1"]
    smiles_out = ["CCO", "c1ccccc1O", "C1CC1OCC"]
    answer = "[OX]"
    st_check_smarts(smiles_in, smiles_out, answer, prompt, idx_offset=3)

    st.markdown(
        """
        ## Ring membership and ring size

        Rings are important in chemistry and many properties differ based on
        whether the atom is in a ring or not. Thus SMARTS allows us to specify
        that an atom is in a ring using the `R` keyword and specify the size of
        the smallest ring an atom is in using the `r` keyword.
        """
    )

    smt = "[CR]=O"
    smiles = ["CC(=O)C", "C1CC(=O)CC1", "c1cccc(C(=O)C2)c12"]
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown(
        """
        As you see, the `R` indicates that the atom must be in a ring. If we

        """
    )

    smt = "[C!R]=O"
    smiles = ["CC(=O)C", "C1CC(=O)CC1", "c1cccc(C(=O)C2)c12"]
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown(
        """
        Here we will touch on another powerful feature of SMARTS, we can specify
        multiple properties at the same time. Here we can combine the `R`, `X`,
        and the `H` keywords to specify that we want cyclic *secondary* amines.
        """
    )

    st.info(
        """
        In SMARTS, `;` is used to mean AND. Here we use `;` to join multiple
        atom properties together.
        """
    )

    smt = "[NR;X3;H1]"
    smiles = ["CNC", "CC=N", "C1N(C)CC1", "C1NCCC1"]
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown(
        """
        The `R` makes sure that the nitrogen is in a ring, the `X3` ensures that
        we only include nitrogens with three neighbors (no imines!), and the
        `H1` ensures that the nitrogen has 1 hydrogen.

        Above, we showed how the `,` keyword can be used to allow flexibility in
        the SMARTS atom matching. However, `,` can also be used to allow
        flexibility in atomic properties as well. For example, we may want
        nitrogens that are in a three- or four-membered ring but not those in a
        five-, six-, or larger-membered ring.
        """
    )

    smt = "[NR;r3,r4]"
    smiles = ["C1NC1", "C1NCC1", "C1NCCC1", "C1NCCCC1"]
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown(
        """
        ### Your turn!
        """
    )

    prompt = "Separate the cyclic ethers from the acyclic ethers."
    smiles_in = ["C1OC1", "C1OC=CC1", "c1cccc2c1OCC2"]
    smiles_out = ["COC", "COC=CC", "c1ccccc1OC"]
    answer = "[OR]"
    st_check_smarts(smiles_in, smiles_out, answer, prompt, idx_offset=4)


def smarts_2_page():
    st.markdown(
        """
        # SMARTS Part 2

        We will continue where we left off in SMARTS Part 1. Previously we
        focused on atomic properties in SMARTS. In this tutorial, we will focus
        on specifying more complex substructures using recursive definitions in
        SMARTS.

        ## Recursion

        If we could only specify atom and bond properties, SMARTS would be
        interesting but somewhat insufficient. Chemical patterns are very
        complex and often, one SMARTS is insufficient to fully capture the
        substructure we are trying to describe. This is where recursion comes
        in. Recursion is like inception, a SMARTS string *within* a SMARTS
        string. First we will go through the basics, then we will see a few
        advanced examples. Then you will really be wowed by the power of SMARTS.

        To start, we'll use a very simple case of recursion which is not very
        useful but is a nice and simple starting point.
        """
    )

    smt = "NC=O"
    smiles = ["CNC(=O)C", "c1ncccc1C(=O)N(C)C", "CNCC", "CNC(=N)C"]
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown(
        """
        As expected, the amides are highlighted while the amine and amidine are
        not. However, with the power of recursion, the SMARTS can be rewritten
        using the `$(...)` format.
        """
    )

    smt = "[N;$(NC=O)]"
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown(
        """
        Again, the amides are selected while the amines and amidines are not.
        However, there is a key difference that only the nitrogen is highlighted
        rather than highlighting the entire amide. Sometimes you do want the
        entire substructure to be selected and sometimes you might only want one
        atom, maybe the reactive atom, to be highlighted.

        Now let's turn it around. Now that we have some power of recursion,
        let's think about how we can use it to specify basic amines. Of course
        writing a SMARTS that can *actually* tell us if an amine is likely basic
        is kind of impossible. However, we can rule out some obvious
        substructures using recursion.
        """
    )

    smiles = ["NC=O", "N=O", "NS(=O)(=O)", "Nc1ncccc1"]
    labels = ["amide", "nitroso", "sulfonamide", "2-aminopyridine"]
    st_display_mol_with_smi(smiles, labels=labels)

    st.markdown(
        """
        We're also primary interested in amines so we'll exclude a few other
        commonly basic nitrogen substructures.
        """
    )

    smiles = ["n1ccccc1", "N=C"]
    labels = ["pyridine", "imine/amidine/guanidine"]
    st_display_mol_with_smi(smiles, labels=labels)

    st.markdown(
        """
        Now we're ready to put it all together. Here is our final SMARTS string.
        """
    )

    exclude = ["NC=[N,O]", "N=O", "NS(=O)(=O)", "Nc1ncccc1", "N=C"]
    smiles = [
        "CNC(=O)C",
        "c1ccccc1[N+](=O)[O-]",
        "NS(=O)(=O)C",
        "CNc1ncccc1",
        "N=CN",
    ]
    smt = f'[N;{";".join([f"!$({smi})" for smi in exclude])}]'
    st_display_mol_with_smi(smiles, substruct=smt)
    smiles = ["CNC", "CN(C)C", "CNC(F)(F)F", "c1ccccc1NC"]
    st_display_mol_with_smi(smiles, substruct=smt, display_substruct=False)

    st.markdown(
        """
        And there we have it, a SMARTS string that's a pretty good filter for
        basic amines. It's probably not good at finding basic amines (poor true
        positive rate) but it's pretty good at excluding non-basic nitrogens (
        good true negative rate).

        This is just the beginning, things can get way crazier! What about
        writing a SMARTS that includes both aldehydes and ketones?
        """
    )

    smt = "[C;$([CH2]=O),$([CH1](=O)[#6]),$([CH0](=O)([#6])[#6])]"
    smiles = ["C=O", "CC=O", "CC(=O)c1ccccc1"]
    st_display_mol_with_smi(smiles, substruct=smt)
    smiles = ["NC(=O)C", "OC(=O)C", "COC(=O)C"]
    st_display_mol_with_smi(smiles, substruct=smt, display_substruct=False)

    st.markdown(
        """
        This SMARTS contains three recursive definitions:
        * [CH2]=O will match formaldehyde, the simplest aldehyde
        * [CH1](=O&rpar;[#6] will match both acetaldehyde and benzaldehyde
        * [CH2](=O&rpar;([#6])[#6] will match all ketones

        Together, these three recursive definitions joined by the `,` keyword
        captures aldehydes and ketones without including carboxylic acids,
        esters, amides, etc.

        ### Your turn!
        """
    )

    prompt = (
        "Separate the esters from the carboxylic acids, amides, carbamates, and"
        " carbonates."
    )
    smiles_in = ["C(=O)OC", "CC(=O)OC", "c1ccccc1C(=O)OC", "CC(=O)Oc1ccccc1"]
    smiles_out = ["CC(=O)O", "CC(=O)NC", "CNC(=O)OC", "COC(=O)OC"]
    answer = "[C;$([CH](=O)O[#6]),$(C([#6])(=O)O[#6])]"
    st_check_smarts(smiles_in, smiles_out, answer, prompt)

    prompt = (
        "Separate the SNAr-reactive aryl fluorides from the non-reactive aryl f"
        "luorides."
    )
    smiles_in = [
        "c1cccnc1F",
        "c1cccc([N+](=O)[O-])c1F",
        "c1cccc(C(=O)OC)c1F",
        "c1cccc(C#N)c1F",
    ]
    smiles_out = ["c1ccccc1F", "c1cccc(N)c1F", "c1ccncc1F", "c1cccc(O)c1F"]
    answer = "[F;$(Fcn),$(Fcc[N+](=O)[O-]),$(FccC(=O)OC),$(FccC#N)]"
    st_check_smarts(smiles_in, smiles_out, answer, prompt, idx_offset=1)
