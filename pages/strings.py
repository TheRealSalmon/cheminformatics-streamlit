import streamlit as st
from PIL import Image

from streamlit_utils import (
    st_display_mol_with_smi,
    st_match_smi_to_mol,
    st_check_smarts,
    st_display_reaction,
    st_run_reaction,
    st_check_smirks,
)

from pages.strings.intro import intro_page
from pages.strings.smiles import smiles_page


def smarts_1_page():
    st.markdown("""
        # SMARTS Part 1

        SMARTS are a SMiles ARbitrary Target Specification. In other words,
        SMARTS are an extension of SMILES that allows even more precise
        specification of chemical structure. While SMILES focuses on
        representing molecules, SMARTS focuses on representing substructures.

        In Part 1, we will focus on introducing the basic syntax which is
        somewhat more complicated. In Part 2, however, we will really get to see
        the power and flexibility of SMARTS.

        ## Why SMARTS? Substructure specification is hard.

        In some cases, it's possible to specify the substructure you have in
        mind with just a SMILES string. For example, quaternary ammoniums can be
        entirely specified with just SMILES.
    """)

    st.info('The highlight indicates a substructure match.')

    smi = 'C[N+](C)(C)C'
    smiles = ['CN', 'CNC', 'CN(C)C', 'C[NH+](C)C']
    st_display_mol_with_smi(smiles, substruct=smi, use_smiles=True)
    smiles = ['C[N+](C)(C)C', 'C[N+](C)(C)C(C)(C)C', 'C[N+](C)(CC1)C1', 'C[N+](C)(C)c1ccccc1']
    st_display_mol_with_smi(smiles, substruct=smi, use_smiles=True, display_substruct=False)

    st.markdown("""
        As it turns out, the SMILES string `C[N+](C)(C)C` is entirely capable of
        differentiating quaternary ammoniums from primary, secondary, and
        tertiary amines.

        However, there are cases where SMILES is insufficient. Consider below
        the case of separating carboxylic acids from esters.
    """)

    smi = 'O=C[OH]'
    smiles = ['O=CO', 'CC(=O)O', 'CC(=O)OC']
    st_display_mol_with_smi(smiles, substruct=smi, use_smiles=True)

    st.markdown("""
        As we see, despite putting an `H` on the oxygen of the carboxylic acid,
        we're unable to separate acids from esters.
    """)

    st.info("""
        This is likely caused by the fact that hydrogens in SMILES are
        second-class citizens in SMILES and are not explicitly matched, even
        when explicitly written in the SMILES. 
    """)

    st.markdown("""        
        We *could* solve this with SMILES by using programming logic to require
        the substructure to match `O=C[OH]` and *not* match `O=COC`. This would
        filter out esters. However, how would we deal with carbamates? What
        about carbonates? Would we have to write a unique filtering logic for each substructure?
    """)
    smiles = ['CNC(=O)OC', 'COC(=O)OC']
    st_display_mol_with_smi(smiles, substruct=smi, use_smiles=True)

    st.markdown("""
        This is where SMARTS come to the rescue. SMARTS are similar to SMILES
        but they are more specific, stricter, and at the same time much more
        flexible.
    """)

    smt = 'O=C[OH]'
    smiles = ['O=CO', 'CC(=O)O', 'CC(=O)OC']
    st_display_mol_with_smi(smiles, substruct=smt)
    smiles = ['CNC(=O)OC', 'COC(=O)OC']
    st_display_mol_with_smi(smiles, substruct=smt, display_substruct=False)

    st.markdown("""
        As we see here, SMARTS is more strict when it comes to counting the
        number of Hs on the oxygen, allowing us to separate carboxylic acids
        from esters.

        Below, we will introduce the SMARTS syntax and hopefully you will be
        convinced that we really do need a *second* molecular string
        representation.

        ## Atoms

        In SMILES we were somewhat restricted in the atomic properties that
        could be specified. However, in SMARTS, the list is much longer:
        * already in SMILES: charge/isotope/# of hydrogens/chirality
        * aromatic/aliphatic
        * degree or # of connections
        * valence (sum of bond orders)
        * if atom is part of a ring
        * the size of the ring the atom is in
        * wildcard

        We won't go through all these subtopics, instead we will focus on:
        * aromatic wildcards, aliphatic wildcards, and general wildcards
        * degree
        * ring membership and ring size

        ## Wildcards

        The most important topic to start with the is the different ways SMARTS
        allows us to represent wildcards. Wildcards allow us to specify some
        flexibility in the molecular structure.

        There are a few types of wildcards in SMARTS. `*` is a general wildcard,
        it will match *any* atom. `A` will match aliphatic atoms while `a` will
        match aromatic atoms. 
    """)

    smt = 'A'
    smiles = ['CNC', 'c1ocnc1', 'C1NCCC1n1cccc1']
    st_display_mol_with_smi(smiles, substruct=smt)

    smt = 'a'
    st_display_mol_with_smi(smiles, substruct=smt)

    smt = '*'
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown("""
        As we see, SMARTS strictly separates aliphatic and aromatic atoms. What
        if we wanted a SMARTS atom that matched both aliphatic and aromatic
        carbons? In that case, we'd need to specify the carbon by using its
        atomic number, `[#6]`.
    """)

    smt = 'C'
    st_display_mol_with_smi(smiles, substruct=smt)

    smt = 'c'
    st_display_mol_with_smi(smiles, substruct=smt)

    smt = '[#6]'
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown("""
        The last type of "wildcard" isn't really a wildcard but is enabled by
        using SMARTS' most powerful feature, logic. We may want a
        "middle-ground" of flexibility in-between specifying atomic number and
        using a wildcard. This is where the `,` character comes in.

        What if we want to match benzenes and pyridines but not anything else?
    """)

    smt = '[c,n]1ccccc1'
    smiles = ['c1ccccc1', 'n1ccccc1', 'n1ncccc1', 'n1cncnc1']
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown("""
        And ta-da, just like that we were able to construct our own "wildcard".
        This is where SMARTS shines, giving us the flexibility to specify really
        whatever any substructure we can think of.

        ### Your turn!
    """)

    prompt = 'Separate the compounds containing nitrogens from the compounds that do not contain nitrogen.'
    smiles_in = ['CNCC', 'c1ncccc1', 'C1CCC1C(=O)NNC']
    smiles_out = ['COCC', 'c1occc1', 'C1CCC1C(=O)OO']
    answer = '[#7]'
    st_check_smarts(smiles_in, smiles_out, answer, prompt)

    prompt = 'Separate the 5-membered aromatic rings from the 6-membered aromatic rings.'
    smiles_in = ['c1occc1', 'c1[nH]ccc1', 'c1sccc1']
    smiles_out = ['c1ccccc1', 'c1ncccc1', 'c1[o+]cccc1']
    answer = 'a1aaaa1'
    st_check_smarts(smiles_in, smiles_out, answer, prompt, idx_offset=1)

    prompt = 'Separate the reactive epoxides and aziridines from the more stable cyclopropanes.'
    smiles_in = ['c1occc1C1OC1', 'C1CCC(N2)C12', 'CC1N(CC(=O)O)C1C']
    smiles_out = ['c1occc1C1CC1', 'C1CCC(C2)C12', 'CC1C(CC(=O)O)C1C']
    answer = 'C1[O,N]C1'
    st_check_smarts(smiles_in, smiles_out, answer, prompt, idx_offset=2)

    st.markdown("""
        Hopefully you feel comfortable working with atoms in SMARTS. However,
        things will only get more and more interesting from here!

        ## Atomic degree

        Let's first start with the example of the quaternary ammonium above. In
        SMILES we can use `C[N+](C)(C)C`. However, this explicitly writes out
        each carbon. With SMARTS, we can specify the degree and the # of
        hydrogens on the `N` to write a more concise substructure string.
    """)

    smt = '[NX4H0+]'
    smiles = ['CN', 'CNC', 'CN(C)C', 'C[NH+](C)C']
    st_display_mol_with_smi(smiles, substruct=smt)
    smiles = ['C[N+](C)(C)C', 'C[N+](C)(C)C(C)(C)C', 'C[N+](C)(CC1)C1', 'C[N+](C)(C)c1ccccc1']
    st_display_mol_with_smi(smiles, substruct=smt, display_substruct=False)

    st.markdown("""
        Breaking it down, the `X4` indicates that this `N` has `4` neighbors.
        The `H0` says there are no hydrogens on this `N` and the `+` gives us
        the charge.

        In another case, `X#` can be used to separate imine-like nitrogens from
        amine-like nitrogens.
    """)

    smt = '[#7X2]'
    smiles = ['CNC', 'c1ccccc1N', 'CC=N', 'c1nc[nH]c1']
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown("""
        ### Your turn!

        Try to solve it using only the `X` keyword.
    """)

    prompt = 'Separate the carbonyls from the ethers/alcohols.'
    smiles_in = ['CC=O', 'C1CNC(=O)C1', 'c1c(=O)[nH]ccc1']
    smiles_out = ['CCO', 'c1ccccc1O', 'C1CC1OCC']
    answer = '[OX1]'
    st_check_smarts(smiles_in, smiles_out, answer, prompt, idx_offset=3)

    st.markdown("""
        ## Ring membership and ring size

        Rings are important in chemistry and many properties differ based on
        whether the atom is in a ring or not. Thus SMARTS allows us to specify
        that an atom is in a ring using the `R` keyword and specify the size of
        the smallest ring an atom is in using the `r` keyword. 
    """)

    smt = '[CR]=O'
    smiles = ['CC(=O)C', 'C1CC(=O)CC1', 'c1cccc(C(=O)C2)c12']
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown("""
        As you see, the `R` indicates that the atom must be in a ring. If we
        wanted the opposite, we could instead use `!`.
    """)

    smt = '[C!R]=O'
    smiles = ['CC(=O)C', 'C1CC(=O)CC1', 'c1cccc(C(=O)C2)c12']
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown("""
        Here we will touch on another powerful feature of SMARTS, we can specify
        multiple properties at the same time. Here we can combine the `R`, `X`,
        and the `H` keywords to specify that we want cyclic *secondary* amines.
    """)

    smt = '[NR;X3;H1]'
    smiles = ['CNC', 'CC=N', 'C1N(C)CC1', 'C1NCCC1']
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown("""
        The `R` makes sure that the nitrogen is in a ring, the `X3` ensures that
        we only include nitrogens with three neighbors (no imines!), and the
        `H1` ensures that the nitrogen has 1 hydrogen.

        Above, we showed how the `,` keyword can be used to allow flexibility in
        the SMARTS atom matching. However, `,` can also be used to allow
        flexibility in atomic properties as well. For example, we may want
        nitrogens that are in a three- or four-membered ring but not those in a
        five-, six-, or larger-membered ring.
    """)

    smt = '[NR;r3,r4]'
    smiles = ['C1NC1', 'C1NCC1', 'C1NCCC1', 'C1NCCCC1']
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown("""
        ### Your turn!
    """)

    prompt = 'Separate the cyclic ethers from the acyclic ethers.'
    smiles_in = ['C1OC1', 'C1OC=CC1', 'c1cccc2c1OCC2']
    smiles_out = ['COC', 'COC=CC', 'c1ccccc1OC']
    answer = '[#6][OR][#6]'
    st_check_smarts(smiles_in, smiles_out, answer, prompt, idx_offset=4)

def smarts_2_page():
    st.markdown("""
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
    """)

    smt = 'NC=O'
    smiles = ['CNC(=O)C', 'c1ncccc1C(=O)N(C)C', 'CNCC', 'CNC(=N)C']
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown("""
        As expected, the amides are highlighted while the amine and amidine are
        not. However, with the power of recursion, the SMARTS can be rewritten
        using the `$(...)` format.
    """)

    smt = '[N;$(NC=O)]'
    st_display_mol_with_smi(smiles, substruct=smt)

    st.markdown("""
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
    """)

    smiles = ['NC=O', 'N=O', 'NS(=O)(=O)', 'Nc1ncccc1']
    labels = ['amide', 'nitroso', 'sulfonamide', '2-aminopyridine']
    st_display_mol_with_smi(smiles, labels=labels)

    st.markdown("""
        We're also primary interested in amines so we'll exclude a few other
        commonly basic nitrogen substructures.
    """)

    smiles = ['n1ccccc1', 'N=C']
    labels = ['pyridine', 'imine/amidine/guanidine']
    st_display_mol_with_smi(smiles, labels=labels)

    st.markdown("""
        Now we're ready to put it all together. Here is our final SMARTS string.
    """)

    exclude = ['NC=[N,O]', 'N=O', 'NS(=O)(=O)', 'Nc1ncccc1', 'N=C']
    smiles = ['CNC(=O)C', 'c1ccccc1[N+](=O)[O-]', 'NS(=O)(=O)C', 'CNc1ncccc1', 'N=CN']
    smt = f'[N;{";".join([f"!$({smi})" for smi in exclude])}]'
    st_display_mol_with_smi(smiles, substruct=smt)
    smiles = ['CNC', 'CN(C)C', 'CNC(F)(F)F', 'c1ccccc1NC']
    st_display_mol_with_smi(smiles, substruct=smt, display_substruct=False)

    st.markdown("""
        And there we have it, a SMARTS string that's a pretty good filter for
        basic amines. It's probably not good at finding basic amines (poor true
        positive rate) but it's pretty good at excluding non-basic nitrogens (
        good true negative rate).

        This is just the beginning, things can get way crazier! What about
        writing a SMARTS that includes both aldehydes and ketones?
    """)

    smt = '[C;$([CH2]=O),$([CH1](=O)[#6]),$([CH0](=O)([#6])[#6])]'
    smiles = ['C=O', 'CC=O', 'CC(=O)c1ccccc1']
    st_display_mol_with_smi(smiles, substruct=smt)
    smiles = ['NC(=O)C', 'OC(=O)C', 'COC(=O)C']
    st_display_mol_with_smi(smiles, substruct=smt, display_substruct=False)

    st.markdown("""
        This SMARTS contains three recursive definitions:
        * [CH2]=O will match formaldehyde, the simplest aldehyde
        * [CH1](=O&rpar;[#6] will match both acetaldehyde and benzaldehyde
        * [CH2](=O&rpar;([#6])[#6] will match all ketones, aliphatic and aromatic

        Together, these three recursive definitions joined by the `,` keyword
        captures aldehydes and ketones without including carboxylic acids,
        esters, amides, etc.

        ### Your turn!
    """)

    prompt = 'Separate the esters from the carboxylic acids, amides, carbamates, and carbonates.'
    smiles_in = ['C(=O)OC', 'CC(=O)OC', 'c1ccccc1C(=O)OC', 'CC(=O)Oc1ccccc1']
    smiles_out = ['CC(=O)O', 'CC(=O)NC', 'CNC(=O)OC', 'COC(=O)OC']
    answer = '[C;$([CH](=O)O[#6]),$(C([#6])(=O)O[#6])]'
    st_check_smarts(smiles_in, smiles_out, answer, prompt)

    prompt = 'Separate the SNAr-reactive aryl fluorides from the non-reactive aryl fluorides.'
    smiles_in = ['c1cccnc1F', 'c1cccc([N+](=O)[O-])c1F', 'c1cccc(C(=O)OC)c1F', 'c1cccc(C#N)c1F']
    smiles_out = ['c1ccccc1F', 'c1cccc(N)c1F', 'c1ccncc1F', 'c1cccc(O)c1F']
    answer = '[F;$(Fcn),$(Fcc[N+](=O)[O-]),$(FccC(=O)OC),$(FccC#N)]'
    st_check_smarts(smiles_in, smiles_out, answer, prompt, idx_offset=1)

def smirks_page():
    st.markdown("""
        # SMIRKS

        As far as I can tell, unlike SMILES or SMARTS, SMIRKS doesn't actually
        stand for anything. That means I get to make something up! How about
        SMarts-Imitating ReaKtion String.` So far we've covered string
        representations of molecules (SMILES) and substructures (SMARTS) but my
        favorite cheminformatic string representation is definitely the SMIRKS.

        ## Atom mapping

        Let's start by looking at one of the simplest reaction in organic
        chemistry, the reaction of methyl bromide with sodium ethoxide.
    """)

    reactants = ['CBr', '[O-]CC']
    products = ['COCC', '[Br-]']
    st_display_reaction(reactants=reactants, products=products)

    st.markdown("""
        In this classic SN2 reaction, a `C-Br` bond is being broken and a new
        `O-C` bond is being formed. As it turns out, SMIRKS has a way to encode
        this very reaction in a string format!

        SMIRKS is very similar to SMARTS. Actually, the only difference between
        SMIRKS and SMARTS are:
        * you must have exactly 2 `>` characters in your string
            * the format is `[reactants]>[reagents]]>[products]
        * you must give the atom mapping (what atom in the reacant is what atom in the product)
            * the format is `:#` where `#` is the atom number
    """)

    smk = '[C:1][Br:2].[O-:3][C:4][C:5]>>[C:1][O-0:3][C:4][C:5].[Br-:2]'
    reactants = ['CBr', '[O-]CC']
    st_run_reaction(smk, reactants)

    st.markdown("""
        Let's break down this SMIRKS:
        * `[C:1][Br:2]` represents methyl bromide
            * the carbon is now atom 1 and the bromine is atom 2
        * `[O-:3][C:4][C:5]` represents the ethoxide
            * the atoms in ethoxide are now assigned an atom number as well
        * `[C:1][O:3][C:4][C:5]` represents the product, ethyl methyl ether
            * as we see in the atom numbers, `[C:1]` comes from the methyl bromide while the rest of the atoms come from ethoxide
        * `[Br-:2]` represents bromide, the by-product of this reaction

        What makes SMIRKS powerful is that it is built off of SMARTS. It's meant
        to be generalizable, and again, capable of expressing our "chemical
        intuition" of what should and shouldn't react.

        For example, we actually don't need to include all the atoms on the
        ethoxide. We don't need both atoms `4` and `5`, we only actually need
        `4`.
    """)

    smk = '[C:1][Br:2].[O-:3][C:4]>>[C:1][O-0:3][C:4].[Br-:2]'
    st_run_reaction(smk, reactants)

    st.markdown("""
        This is because atom `5` is not involved in the reaction. We should also
        rewrite this SMIRKS to ignore the by-product bromide ion.
    """)

    smk = '[C:1]Br.[O-:3][C:4]>>[C:1][O-0:3][C:4]'
    st_run_reaction(smk, reactants)

    st.markdown("""
        And there we have it, a SMIRKS string that does a pretty good job at
        encoding a basic SN2 reaction.

        ### Your turn!

        I know that writing SMIRKS is somewhat more tedious than SMILES or
        SMARTS. But I would strongly encourage you to give these a serious try
        before reaching for the answer button.
    """)

    prompt = 'Reaction: Etherification'
    reactantss = [
        ['CCl', 'Oc1ccccc1']
    ]
    productss = [
        ['COc1ccccc1']
    ]
    answer = '[C:1]Cl.[O:2][c:3]>>[C:1][O:2][c:3]'
    st_check_smirks(reactantss, productss, answer, prompt)

    prmopt = 'Reaction: Reductive amination'
    reactantss = [
        ['CC=O', 'NC1CC1']
    ]
    productss = [
        ['CCNC1CC1']
    ]
    answer = '[C:1]=O.[N:2][C:3]>>[C:1][N:2][C:3]'
    st_check_smirks(reactantss, productss, answer, prompt, idx_offset=1)

    st.markdown("""
        ## Making our SMIRKS more smart

        What makes SMIRKS really powerful is the power of SMARTS that it
        inherits. For example, what if we wanted a SMIRKS that could capture
        reaction of alkoxides with methyl chloride and bromide?
    """)

    smk = '[C:1][Cl,Br].[O:2][C:3]>>[C:1][O:2][C:3]'
    reactants = ['CCl', 'OC']
    st_run_reaction(smk, reactants)
    st.markdown('---')
    reactants = ['CBr', 'OC']
    st_run_reaction(smk, reactants, display_smk=False)

    st.markdown("""
        As we see, the key component here is the `[Cl,Br]`. This gives the
        SMIRKS some flexibility to select either methyl chloride or bromide.

        What if we wanted to go the other way and add some restrictions to our
        SMIRKS? For example this SMIRKS allows tert-butyl alcohol to react but
        we know that SN2 is less likely for hindered nucleophiles.
    """)

    smk = '[C:1][Cl,Br].[O:2][C:3]>>[C:1][O:2][C:3]'
    reactants = ['CCBr', 'OC']
    st_run_reaction(smk, reactants)
    st.markdown('---')
    reactants = ['CCBr', 'OC(C)(C)C']
    st_run_reaction(smk, reactants, display_smk=False)

    st.markdown("""
        What if we wanted to exclude t-butyl alcohols? We could add in some
        SMARTS properties to `[C:3]`.
    """)

    smk = '[C:1][Cl,Br].[O:2][C;!H0:3]>>[C:1][O:2][C:3]'
    reactants = ['CCBr', 'OC']
    st_run_reaction(smk, reactants)
    st.markdown('---')
    reactants = ['CCBr', 'OC(C)(C)C']
    st_run_reaction(smk, reactants, display_smk=False)

    st.markdown("""
        Specifying the carbon attached to the alcohol as `[C;!H0:3]` excludes
        the t-butyl alcohol as that partcular carbon has no hydrogens. However,
        it's not that simple, what about phenols? Our SMIRKS excludes phenols
        but we know that phenols should be able to react.
    """)

    reactants = ['CCBr', 'Oc1ccccc1']
    st_run_reaction(smk, reactants)

    st.markdown("""
        There are two problems.
        * We are specifying the carbon attached to the oxygen as `C`. This will exclude aromatic carbons, like the carbons in the phenol ring.
        * We are requiring the carbon attached to the oxygen to have 1, 2, or 3 hydrogens but the carbon attached to the phenol oxygen has 0.

        This can be fixed by using `#6` in place of `C` and using recursion to
        aliphatic carbon logic from aromatic carbon logic.
    """)

    smk = '[C:1][Cl,Br].[O:2][#6;$([C!H0]),$(c):3]>>[C:1][O:2][#6:3]'
    reactants = ['CCBr', 'OC']
    st_run_reaction(smk, reactants)
    st.markdown('---')
    reactants = ['CCBr', 'OC(C)(C)C']
    st_run_reaction(smk, reactants, display_smk=False)
    st.markdown('---')
    reactants = ['CCBr', 'Oc1ccccc1']
    st_run_reaction(smk, reactants, display_smk=False)

    st.markdown("""
        Hopefully you see now why SMIRKS is my personal favorite. Using the
        flexibility and specificity of SMARTS, we're able to encode the logic of
        a chemical reaction into a sequence of letters, numbers, and symbols.

        ### Your turn!
    """)

    prompt = 'Reaction: SNAr of aryl fluorides.'
    reactantss = [
        ['c1cccnc1F', 'N1CCCC1'],
        ['c1ccccc1F', 'N1CCCC1'],
    ]
    productss = [
        ['c1cccnc1N1CCCC1'],
        ['']
    ]
    answer = 'F[c:1][n:2].[N:3]>>[N:3][c:1][n:2]'
    st_check_smirks(reactantss, productss, answer, prompt, idx_offset=2)

    prompt = 'Reaction: Boc deprotection.'
    reactantss = [
        ['C1CCN1C(=O)OC(C)(C)C'],
    ]
    productss = [
        ['C1CCN1'],
    ]
    answer = '[#7:1]C(=O)OC(C)(C)C>>[#7:1]'
    st_check_smirks(reactantss, productss, answer, prompt, idx_offset=3)

    prompt = 'Reaction: Reduction of ketones.'
    reactantss = [
        ['C1C(=O)CCC1'],
        ['C1C(=O)NCC1'],
        ['C1C(=O)OCC1'],
    ]
    productss = [
        ['C1C(O)CCC1'],
        [''],
        [''],
    ]
    answer = '[C;$(C(C)C):1]=[O:2]>>[C:1][O:2]'
    st_check_smirks(reactantss, productss, answer, prompt, idx_offset=4)

    st.markdown("""
        And there we have it, a primer on SMIRKS. Once you have enough practice
        with SMIRKS, many doors will open. SMIRKS is used for many, many tasks
        in cheminformatics:
        * exchanging between different tautomers
        * virtually deprotecting molecules, useful for trawling commercial building blocks
        * building virtual compound libraries
    """)


subpages = [
    ('0 - intro', intro_page),
    ('1 - SMILES', smiles_page),
    ('2 - SMARTS Part 1', smarts_1_page),
    ('3 - SMARTS Part 2', smarts_2_page),
    ('4 - SMIRKS Part 1', smirks_page),
]
st.sidebar.markdown('## sub-topics')
subpage = st.sidebar.selectbox(
    'select page',
    subpages,
    format_func=lambda x: x[0]
)
st.sidebar.markdown("---")
subpage[1]()
