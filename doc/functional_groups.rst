=================
Functional Groups
=================

FGUtils provides a class :py:class:`~fgutils.query.FGQuery` to query a
molecules functional groups. 

Functional group tree
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
