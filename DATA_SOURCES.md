# Data sources and provenance

SOMatic operates on atom-mapped substrate–enzyme records. This repository does not redistribute source-derived labels or structures. Reproduction should preserve each source's original identifiers, stereochemistry, atom mapping, candidate mask, positive-site labels, enzyme family/subtype, and citation.

## Publicly described evaluation sources

| Name used in SOMatic | Primary source | Persistent identifier | Role |
|---|---|---|---|
| EBoMD2-derived / CyProduct | Tian et al., *Journal of Chemical Information and Modeling* (2021) | DOI: [10.1021/acs.jcim.1c00144](https://doi.org/10.1021/acs.jcim.1c00144) | CYP site-of-metabolism structures and labels |
| AZ120 | Chen et al., *Molecular Pharmaceutics* (2025) | DOI: [10.1021/acs.molpharmaceut.5c00740](https://doi.org/10.1021/acs.molpharmaceut.5c00740) | Exact-site and regional external endpoints; keep endpoints separate |
| PublicSOM (exploratory assembly) | GLORY reference: de Bruyn Kops et al., *Frontiers in Chemistry* (2019) | DOI: [10.3389/fchem.2019.00402](https://doi.org/10.3389/fchem.2019.00402) | Exploratory public-source panel; not an independent benchmark |

The corrected working collection contains CYP, UGT, FMO, and AOX1 records assembled from source publications and databases. Some source records cannot be redistributed under an open repository license. Obtain them from the cited primary sources and respect their terms.

## Representation source

- Uni-Mol2-84M: Ji et al., “Exploring molecular pretraining model at scale,” NeurIPS 2024, DOI [10.52202/079017-1489](https://doi.org/10.52202/079017-1489).
- Feature software used by SOMatic: `unimol_tools` 0.1.4.post1.
- Gasteiger charges: Gasteiger and Marsili, *Tetrahedron* (1980), DOI [10.1016/0040-4020(80)80168-2](https://doi.org/10.1016/0040-4020(80)80168-2).

## Required record manifest

For a numerically auditable reproduction, publish or obtain a versioned manifest with these columns:

- `record_id`, `entity_id`, `source_id`;
- canonical isomeric SMILES and the ordered molecular graph;
- `enzyme_family` and original enzyme label;
- zero-based `candidate_atoms` and `positive_atoms`;
- label-completeness and mapping-status fields;
- substrate/scaffold split and association-group identifiers;
- SHA256 checksums for the record table, split manifest, and frozen atom-feature cache.

All heavy atoms are candidates for CYP. For UGT, FMO, and AOX1, candidate sets follow source-listed potential sites; hydrogen labels are mapped to their attached heavy atom. Do not infer missing non-CYP candidate lists from products without documenting that change.

## Reproduction boundary

The public code fixes the model definition, seeds, weighting, feature dimensions, and software environment. Exact reproduction of historical scores additionally depends on the corrected record/split manifest and the frozen Uni-Mol2 feature cache. Historical coordinates and fallback flags were unavailable for some retained cache rows, so a fresh encoder run can reproduce the declared pipeline but is not guaranteed to be bitwise identical to every historical feature row.
