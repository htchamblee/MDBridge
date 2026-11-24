# Standard residue names

`standard_residues.tsv` lists the residue names the generator treats as standard. Facts taken from:

- **Amino acids (20):** IUPAC-IUB Joint Commission on Biochemical Nomenclature, *Nomenclature and
  Symbolism for Amino Acids and Peptides* (1983), 3AA-1 Table 1, three-letter symbols upper-cased
  as in PDB files. https://iupac.qmul.ac.uk/AminoAcid/AA1n2.html, retrieved 2026-11-24.
- **Sec / Pyl:** selenocysteine (IUPAC-IUB JCBN Newsletter 1999, symbol Sec/U) and pyrrolysine
  (Pyl/O), with their wwPDB CCD ids SEC and PYL (CC0).
- **Nucleotides (10):** wwPDB *Atomic Coordinate Entry Format Version 3.3*, section 3 (SEQRES):
  ribo- forms A, C, G, U, I and deoxy- forms DA, DC, DG, DT, DI.
  https://www.wwpdb.org/documentation/file-format-content/format33/sect3.html, retrieved 2026-11-24.
- **Water:** HOH, as in the HETATM examples of the same document, section 9 (and CCD id HOH).

**Licence:** residue names and symbols are facts, not copyrightable expression; the IUPAC
recommendations and the wwPDB format guide (© wwPDB) are cited, not copied. The atom-column rule
used by the PDB writer ("Alignment of one-letter atom name such as C starts at column 14, while
two-letter atom name such as FE starts at column 13", ATOM record, section 9) is from the same guide.
