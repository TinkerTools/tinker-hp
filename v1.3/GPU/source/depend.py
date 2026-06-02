#!/usr/bin/env python3

import sys
import re
from pathlib import Path

# ==========================================================
# arguments
# ==========================================================

if len(sys.argv) != 2 or sys.argv[1] not in ["cpu", "host", "gpu", "device"]:
    print("usage: depend.py [host|cpu|device|gpu]")
    sys.exit(1)

mode = "h" if sys.argv[1] in ["cpu","host"] else "d"

# ==========================================================
# source directory
# ==========================================================

source_dir = Path(__file__).resolve().parent

# ==========================================================
# structures
# ==========================================================

dep_dict = {}      # clé logique -> objet ou fichier
obj_source = {}    # objet -> fichier source
obj_deps = {}      # objet -> dépendances

# ==========================================================
# regex robustes
# ==========================================================

re_module = re.compile(
    r'^\s*module\s+(?!procedure)(\w+)',
    re.IGNORECASE
)

re_use = re.compile(
    r'^\s*use\s*(?:,\s*intrinsic\s*::)?\s*(\w+)',
    re.IGNORECASE
)

re_include = re.compile(
    r'#include\s*[<"]([^">]+)[">]'
)

# ==========================================================
# classification fichiers
# ==========================================================

def classify_file(file):

    name = file.name

    # inclusion fortran
    if name.endswith(".inc.f90") or name.endswith(".tpl.f90") or name.endswith(".h.f90"):
        return None

    # header C++
    if name.endswith(".h"):
        return None

    # CUDA Fortran
    if name.endswith(".cu.f90"):
        if mode == "h":
            return "ignore"
        base = name[:-7]
        return f"cuf.{base}.o"

    # module Fortran
    if name.startswith("MOD_") and name.endswith(".f90"):
        return name.replace(".f90", ".o")

    # Fortran classique
    if name.endswith(".f90"):
        base = name[:-4]
        return f"{base}.o"

    # CUDA C++
    if name.endswith(".cu"):
        if mode == "h":
            return "ignore"
        base = name[:-3]
        return f"cu.{base}.o"

    # C++
    if name.endswith(".cpp"):
        if mode == "h":
            return "ignore"
        base = name[:-4]
        return f"cc.{base}.o"

    return "ignore"

# ==========================================================
# ETAPE 1 : scan fichiers
# ==========================================================

for file in source_dir.iterdir():

    if not file.is_file():
        continue

    obj = classify_file(file)

    if obj == "ignore":
        continue

    # fichiers d'inclusion
    if obj is None:
        dep_dict[file.name] = str(source_dir / file.name)
        continue

    # modules Fortran
    if (file.name.startswith("MOD_") and file.suffix == ".f90")\
     or file.name.endswith(".cu.f90"):

        module_name = None

        with open(file, errors="ignore") as f:
            for line in f:
                m = re_module.match(line)
                if m:
                    module_name = m.group(1).lower()
                    break

        if module_name:
            dep_dict[module_name] = obj
        else:
            dep_dict[file.name] = obj

    else:
        dep_dict[file.name] = obj

    obj_source[obj] = file
    obj_deps[obj] = set()


#  Added Extra dependencies
#  ( Following Files containing multiple modules )
#  ------------------------------------------------------------------------
dep_dict["tinmemory"]  = "MOD_memory.o"
dep_dict["potent"]     = "MOD_potent.o"
dep_dict["polar_temp"] = "MOD_polar.o"


# ==========================================================
# ETAPE 2 : dépendances
# ==========================================================

for obj, file in obj_source.items():

    with open(file, errors="ignore") as f:

        for line in f:

            # USE modules
            m = re_use.match(line)
            if m:
                mod = m.group(1).lower()
                if mod in dep_dict and dep_dict[mod] != obj:
                    obj_deps[obj].add(dep_dict[mod])

            # includes
            for inc in re_include.findall(line):

                inc = Path(inc).name

                if inc in dep_dict and dep_dict[inc] != obj:
                    obj_deps[obj].add(dep_dict[inc])

# ==========================================================
# ETAPE 3 : écriture Makefile
# ==========================================================

with open("./depend.mk", "w") as out:

    out.write("# automatically generated\n\n")

    for obj in sorted(obj_source):

        deps = sorted(obj_deps[obj])

        if not deps:
            continue

        dep_line = " ".join(deps)

        out.write(f"{obj}: {dep_line}\n")
