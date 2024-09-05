=================
Functional Groups
=================

FGUtils provides a class :py:class:`~fgutils.query.FGQuery` to retrieve a
molecules functional groups. It can be used directly with the preconfigured
list of functional groups as listed in `Functional Group Tree`_ or by
specifying your own patterns using :py:class:`~fgutils.fgconfig.FGConfig`
and :py:class:`~fgutils.fgconfig.FGConfigProvider`.

Get functional groups in molecule
=================================

Common functional groups of a molecule can be retrieved by the
:py:class:`~fgutils.query.FGQuery` class. The query can directly use the
molecules SMILES or the molecular graph representation of the compound as
networkx graph. The following example demonstrates how to get the functional
groups from *acetylsalicylic acid*::

    >>> from fgutils import FGQuery

    >>> smiles = "O=C(C)Oc1ccccc1C(=O)O"  # acetylsalicylic acid
    >>> query = FGQuery()
    >>> query.get(smiles)
    [("ester", [0, 1, 3]), ("carboxylic_acid", [10, 11, 12])]


Get changing groups in reaction
===============================

The extended :ref:`graph-syntax` enables the description of reaction
mechanisms by specifying bond changes in ``<>`` brackets. Functional group
patterns can therefor also specify bond changes. Querying bond changes can be
used to look for a changing functional groups in a reaction. The following
example demonstrates how to check for a nucleophilic addition-elimination
reaction on a carbonyl group::

    >>> from fgutils.query import FGQuery, FGConfig

    >>> smiles = "[C:1][C:2](=[O:3])[O:4][C:5].[O:6]>>[C:1][C:2](=[O:3])[O:6].[O:4][C:5]"
    >>> fgconfig = FGConfig(name="carbonyl-AE", pattern="C(=O)(<0,1>R)<1,0>R")
    >>> query = FGQuery(config=fgconfig, require_implicit_hydrogen=False)
    >>> query.get(smiles)
    [("carbonyl-AE", [2, 3, 4, 6])]


Functional Group Tree
=====================

.. code-block::

    Functional Group                    Parents                  Pattern
    --------------------------------------------------------------------------
    ether                               [ROOT]                   ROR
    ├── ketal                           [ether]                  RC(OR)(OR)R
    │   ├── acetal                      [ketal]                  RC(OC)(OC)H
    │   └── hemiketal                   [ketal, alcohol]         RC(OH)(OR)R
    │       └── hemiacetal              [hemiketal]              RC(OC)(OH)H
    ├── ester                           [ketone, ether]          RC(=O)OR
    │   ├── anhydride                   [ester]                  RC(=O)OC(=O)R
    │   ├── peroxy_acid                 [ester, peroxide]        RC(=O)OOH
    │   ├── carbamate                   [ester, amide]           ROC(=O)N(R)R
    │   └── carboxylic_acid             [ester, alcohol]         RC(=O)OH
    ├── alcohol                         [ether]                  COH
    │   ├── hemiketal                   [ketal, alcohol]         RC(OH)(OR)R
    │   │   └── hemiacetal              [hemiketal]              RC(OC)(OH)H
    │   ├── carboxylic_acid             [ester, alcohol]         RC(=O)OH
    │   ├── enol                        [alcohol]                C=COH
    │   ├── primary_alcohol             [alcohol]                CCOH
    │   │   └── secondary_alcohol       [primary_alcohol]        C(C)(C)OH
    │   │       └── tertiary_alcohol    [secondary_alcohol]      C(C)(C)(C)OH
    │   └── phenol                      [alcohol]                C:COH
    └── peroxide                        [ether]                  ROOR
        └── peroxy_acid                 [ester, peroxide]        RC(=O)OOH
    thioether                           [ROOT]                   RSR
    └── thioester                       [ketone, thioether]      RC(=O)SR
    amine                               [ROOT]                   RN(R)R
    ├── amide                           [ketone, amine]          RC(=O)N(R)R
    │   └── carbamate                   [ester, amide]           ROC(=O)N(R)R
    └── anilin                          [amine]                  C:CN(R)R
    carbonyl                            [ROOT]                   C(=O)
    ├── ketene                          [carbonyl]               RC(R)=C=O
    └── ketone                          [carbonyl]               RC(=O)R
        ├── amide                       [ketone, amine]          RC(=O)N(R)R
        │   └── carbamate               [ester, amide]           ROC(=O)N(R)R
        ├── thioester                   [ketone, thioether]      RC(=O)SR
        ├── ester                       [ketone, ether]          RC(=O)OR
        │   ├── anhydride               [ester]                  RC(=O)OC(=O)R
        │   ├── peroxy_acid             [ester, peroxide]        RC(=O)OOH
        │   ├── carbamate               [ester, amide]           ROC(=O)N(R)R
        │   └── carboxylic_acid         [ester, alcohol]         RC(=O)OH
        ├── acyl_chloride               [ketone]                 RC(=O)Cl
        └── aldehyde                    [ketone]                 RC(=O)H
    nitrose                             [ROOT]                   RN=O
    └── nitro                           [nitrose]                RN(=O)O
    nitrile                             [ROOT]                   RC#N
