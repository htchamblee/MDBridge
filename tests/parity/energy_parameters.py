"""The energy analysis' parameters without parm7 (docs/DIFFERENCES.md AN-25 - AN-27, CIT-04), for
MDATK's patched energy analysis (improvements_analyses.py): Lennard-Jones parameters from Amber
force-field residue templates, else from Amber atom types, and partial charges from completely
matched templates. The same rules as src/analysis/energy_parameters.cpp (MDBridge's own design,
not MDAnalysis code), on the tables tools/gen_energy_lj.py writes for it (load()).

energy_parameters(u, indices, atom_elements, guess_bonds) -> dict with per-atom lists rstar,
epsilon, charge, layer ("template" / "typed"), template_charge, and force_fields, references,
default_elements, default_atoms, chargeless_residues, as EnergyParameters in C++.
"""

import importlib.util
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

HYDROGEN_REACH = 1.45
MOLECULE_TYPES = {"C": "C", "N": "N", "O": "O", "P": "P", "S": "S", "F": "F", "CL": "Cl", "BR": "Br", "I": "I"}
ION_TYPES = {"LI": "Li+", "NA": "Na+", "K": "K+", "RB": "Rb+", "CS": "Cs+", "F": "F-", "CL": "Cl-", "BR": "Br-",
             "I": "I-", "MG": "Mg2+", "CA": "Ca2+", "SR": "Sr2+", "BA": "Ba2+", "MN": "Mn2+", "FE": "Fe2+",
             "CO": "Co2+", "NI": "Ni2+", "CU": "Cu2+", "ZN": "Zn2+", "CD": "Cd2+"}

_TABLES = None


def tables():
    """tools/gen_energy_lj.py's load() (the tree's own tools/, as regenerate_all.py copies it)."""
    global _TABLES
    if _TABLES is None:
        spec = importlib.util.spec_from_file_location("gen_energy_lj", os.path.join(ROOT, "tools", "gen_energy_lj.py"))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        data = module.load()
        data["by_name"] = {t["name"]: k for k, t in enumerate(data["templates"])}
        _TABLES = data
    return _TABLES


def short_citation(key):
    return tables()["short"].get(key, key)


def _strip_charge(text):
    return text.rstrip("+-0123456789")


def _has_variants(data, t):
    return any(o["variant"] and o["base"] == t["name"] and o["forcefield"] == t["forcefield"]
               for o in data["templates"])


def _add_with_variants(data, k, out):
    """The template, and its terminal variants when it is the middle form; returns filterable."""
    out.append(k)
    t = data["templates"][k]
    if t["variant"] or not _has_variants(data, t):
        return False
    for v, o in enumerate(data["templates"]):
        if o["variant"] and o["base"] == t["name"] and o["forcefield"] == t["forcefield"]:
            out.append(v)
    return True


def _map_atoms(data, t, names, elements):
    """(map, complete) or None: as mapAtoms in energy_parameters.cpp."""
    atoms = t["atoms"]
    index = {}
    for a, atom in enumerate(atoms):
        index.setdefault(atom[0], a)
    mapping = [-1] * len(names)
    used = [False] * len(atoms)
    for k, name in enumerate(names):
        a = index.get(name, -1)
        if a >= 0:
            if used[a]:
                return None
            used[a] = True
            mapping[k] = a
    aliases = data["atom_aliases"]
    for k, name in enumerate(names):
        if mapping[k] >= 0:
            continue
        target = aliases.get(t["specific"], {}).get(name) if t["specific"] else None
        if target is None and t["generic"]:
            target = aliases.get(t["generic"], {}).get(name)
        a = index.get(target, -1) if target is not None else -1
        if a < 0 or used[a]:
            return None
        used[a] = True
        mapping[k] = a
    for k in range(len(names)):
        element = atoms[mapping[k]][1]
        if elements[k] and element and elements[k] != element:
            return None
    return mapping, all(used)


def energy_parameters(u, indices, atom_elements, guess_bonds):
    """indices: atom indices (A then B). atom_elements(AtomGroup) -> upper-case elements (MDATK's
    patched selections); guess_bonds(AtomGroup) -> local index pairs (AN-14), on the current frame."""
    data = tables()
    templates = data["templates"]
    n_atoms = u.atoms.n_atoms
    elements = [str(e) for e in atom_elements(u.atoms)]
    names = [str(x) for x in u.atoms.names]
    resindex = [int(r) for r in u.atoms.resindices]
    n_res = len(u.residues)
    resnames = [str(x) for x in u.residues.resnames]
    res_atoms = [[] for _ in range(n_res)]
    for a, r in enumerate(resindex):
        res_atoms[r].append(a)
    positions = u.atoms.positions  # float32, the current frame
    try:
        bonds = [(int(x), int(y)) for x, y in u.bonds.indices]
    except Exception:  # noqa: BLE001 - no bond information
        bonds = []
    partners = [[] for _ in range(n_atoms)]
    res_bonded = [False] * n_res
    for x, y in bonds:
        partners[x].append(y)
        partners[y].append(x)
        res_bonded[resindex[x]] = True
        res_bonded[resindex[y]] = True
    guessed = {"done": False, "partners": None}
    symbols = set(data["element_symbols"])
    vsites = set(data["virtual_sites"])
    vnames = {n for _, n in vsites}

    def is_extra_point(a):
        name = names[a].strip().upper()
        if (resnames[resindex[a]].strip().upper(), name) in vsites:
            return True
        if elements[a] in symbols:
            return False
        return name == "MW" or name.startswith("EP") or name.startswith("LP") or name in vnames

    extra = [is_extra_point(a) for a in range(n_atoms)]

    def guessed_partners(atom):
        if not guessed["done"]:
            guessed["done"] = True
            keep = [False] * n_res
            for a in indices:
                r = resindex[a]
                for q in (r - 1, r, r + 1):
                    if 0 <= q < n_res:
                        keep[q] = True
            group = [a for a in range(n_atoms) if keep[resindex[a]]]
            gp = [[] for _ in range(n_atoms)]
            for i, j in guess_bonds(u.atoms[group]):
                i, j = group[int(i)], group[int(j)]
                gp[i].append(j)
                gp[j].append(i)
            guessed["partners"] = gp
        return guessed["partners"][atom]

    def linked(r1, r2):
        if r1 < 0 or r2 < 0 or r1 >= n_res or r2 >= n_res:
            return False
        topo = res_bonded[r1] or res_bonded[r2]
        for a in res_atoms[r1]:
            for b in (partners[a] if topo else guessed_partners(a)):
                if resindex[b] == r2:
                    return True
        return False

    def neighbours(atom):
        out = partners[atom] if res_bonded[resindex[atom]] else guessed_partners(atom)
        return [b for b in out if not extra[b]]

    atom_template = [-1] * n_atoms
    atom_template_atom = [-1] * n_atoms
    atom_charged = [False] * n_atoms

    def assign(atoms, k, mapping, charged):
        for i, a in enumerate(atoms):
            atom_template[a] = k
            atom_template_atom[a] = mapping[i]
            atom_charged[a] = charged and len(atoms) == len(res_atoms[resindex[a]])

    def match_ion(r, element):
        atom = res_atoms[r][0]
        resname = resnames[r].strip()
        name = names[atom].strip()
        for key in (resname, resname.upper(), _strip_charge(resname).upper(), name, name.upper(),
                    _strip_charge(name).upper()):
            if not key:
                continue
            k = data["by_name"].get(key)
            if k is not None and templates[k]["family"] == "ion" and templates[k]["atoms"][0][1] == element:
                return k
        return None

    def match_residue(r):
        atoms = [a for a in res_atoms[r] if not extra[a]]
        extra_points = len(atoms) < len(res_atoms[r])
        if not atoms:
            return
        rnames = [names[a] for a in atoms]
        relements = [elements[a] for a in atoms]
        resname = resnames[r].strip().upper()
        if len(atoms) == 1 and not extra_points:
            ion = match_ion(r, relements[0])
            if ion is not None:
                assign(atoms, ion, [0], True)
                return
        candidates = []
        filterable = False
        k = data["by_name"].get(resname)
        if k is not None and templates[k]["family"] != "ion":
            filterable = _add_with_variants(data, k, candidates) or filterable
        else:
            for target in data["residue_aliases"].get(resname, []):
                kt = data["by_name"].get(target)
                if kt is not None:
                    filterable = _add_with_variants(data, kt, candidates) or filterable
        if not candidates:
            return
        containing = []
        for k in candidates:
            m = _map_atoms(data, templates[k], rnames, relements)
            if m is not None:
                containing.append((k, m[0], m[1]))
        completes = [c for c in containing if c[2]]
        if len(completes) == 1:
            assign(atoms, completes[0][0], completes[0][1], True)
            return
        if filterable:
            start = not linked(r - 1, r)
            end = not linked(r, r + 1)

            def keep(c):
                t = templates[c[0]]
                if not t["variant"] and not _has_variants(data, t):
                    return True
                if t["family"] == "protein":
                    wanted = "-" if start and end else "N" if start else "C" if end else ""
                else:
                    wanted = "N" if start and end else "5" if start else "3" if end else ""
                return t["variant"] == wanted

            containing = [c for c in containing if keep(c)]
        if not containing:
            return
        first = containing[0]
        for c in containing:
            for i in range(len(atoms)):
                a = templates[c[0]]["atoms"][c[1][i]]
                b = templates[first[0]]["atoms"][first[1][i]]
                if a[4] != b[4] or a[5] != b[5]:
                    return
        assign(atoms, first[0], first[1], False)

    done = [False] * n_res
    for a in indices:
        r = resindex[a]
        if not done[r]:
            done[r] = True
            match_residue(r)

    def hydrogen_parent(h):
        bonded = [b for b in partners[h] if not extra[b]]
        if bonded:
            return min(bonded)
        r = resindex[h]
        best, best_d = -1, 0.0
        for a in res_atoms[r]:
            if a == h or elements[a] == "H" or extra[a]:
                continue
            dx = float(positions[a][0]) - float(positions[h][0])
            dy = float(positions[a][1]) - float(positions[h][1])
            dz = float(positions[a][2]) - float(positions[h][2])
            d = math.sqrt((dx * dx + dy * dy) + dz * dz)
            if d <= HYDROGEN_REACH and (best < 0 or d < best_d):
                best, best_d = a, d
        return best

    out = {"rstar": [], "epsilon": [], "charge": [], "layer": [], "template_charge": [], "force_fields": [],
           "references": [], "default_elements": [], "default_atoms": 0, "chargeless_residues": 0}

    def add(lst, item):
        if item and item not in lst:
            lst.append(item)

    def typed(atom):
        element = elements[atom]
        if extra[atom]:
            return "EP"
        if element == "H":
            parent = hydrogen_parent(atom)
            pe = elements[parent] if parent >= 0 else "N"
            return "H" if pe == "N" else "HO" if pe == "O" else "HS" if pe == "S" else "H1"
        r = resindex[atom]
        if len(res_atoms[r]) == 1 and element in ION_TYPES:
            return ION_TYPES[element]
        if element in MOLECULE_TYPES:
            if element == "C" and len(neighbours(atom)) == 4:
                return "CT"
            if element == "O" and len(neighbours(atom)) == 2:
                hydrogens = sum(1 for b in neighbours(atom) if elements[b] == "H")
                return "OW" if hydrogens == 2 else "OH" if hydrogens == 1 else "OS"
            return MOLECULE_TYPES[element]
        if element in ION_TYPES:
            return ION_TYPES[element]
        add(out["default_elements"], element)
        out["default_atoms"] += 1
        return "C"

    chargeless = set()
    for a in indices:
        if atom_template[a] >= 0:
            t = templates[atom_template[a]]
            ta = t["atoms"][atom_template_atom[a]]
            out["rstar"].append(ta[4])
            out["epsilon"].append(ta[5])
            out["layer"].append("template")
            if atom_charged[a]:
                out["charge"].append(ta[3])
                out["template_charge"].append(True)
            else:
                out["charge"].append(0.0)
                out["template_charge"].append(False)
            add(out["force_fields"], t["forcefield"])
            for key in data["reference_keys"]:  # the C++ table's mask order
                if key in t["references"]:
                    add(out["references"], key)
        else:
            rstar, epsilon, ref = data["types"][typed(a)]
            out["rstar"].append(rstar)
            out["epsilon"].append(epsilon)
            out["layer"].append("typed")
            out["charge"].append(0.0)
            out["template_charge"].append(False)
            add(out["references"], ref)
        if not out["template_charge"][-1]:
            chargeless.add(resindex[a])
    out["chargeless_residues"] = len(chargeless)
    return out
