#
#     Sorbonne University
#     Washington University in Saint Louis
#     University of Texas at Austin
#     
#     Manage Fortran Dependencies for Tinker-HP
#
#     This script produces Makefile rules to be included in the Makefile 
#     when configure generates it
#     
#     Luc-Henri Jolly (2019)

import re, sys, os, getopt

#
# gets the dependencies for each FORTRAN source 
#
def get_deps_and_mods(filename):
    deps = []
    mods = []
    f = open(filename)
    if not f:
       print ("ERROR: unable to open %s%s" % filename)
       sys.exit(1)
#
# Match  USE statement in fortran sources (case insensitive)
#
    use_line_re  = re.compile("(?i)^\s*use\s+(\S.+?)\s*(?:,\s*only\s*:.*)?$")
#
# Match  other modules names in the same USE statement on a continuation line
# (i.e.  use modulea,
#       &    moduleb
#
    cont_line_re = re.compile("^(.*)&\s*$")
#
# Match module statement in fortran sources (case insensitive)
#
    mod_line_re  = re.compile("(?i)^\s*module\s+(\S+)\s*$")
#
# Split USE or MODULE statement lines in pieces
#
    split_re     = re.compile("\s*,\s*")
#
# Match any MODULE or USE name
# (except MPI and iso_c_binding, which are not in a source file,  but come with the compiler)
#
    dep_re       = re.compile("(?i)(^(?!.*mpi|.*iso_c_binding).*$)")
    mod_re       = re.compile("(?i)(^(?!.*mpi|.*iso_c_binding).*$)")
    within_use_statement = False
    for line in f:
        match = use_line_re.search(line)
        if match:
           within_use_statement = True
           rest_line = match.group(1)
        else:
           rest_line = line
        if within_use_statement:
#
# We're in USE statement
#
           match = cont_line_re.search(rest_line)
           if match:
               rest_line = match.group(1)
           else:
               within_use_statement = False
               line_items = split_re.split(rest_line.strip())
#
# Loop over all modules found in USE statements
#
           for item in line_items:
               if item:
                  match = dep_re.match(item)
                  if match:
#
# Add the MOD_file.o corresponding to the modules found. Just one time
#
                     dep = match.expand("MOD_\\1.o")
                     if dep not in deps:
                        deps.append(dep)
        else:
#
# We're in a MODULE statement
#
           match = mod_line_re.search(line)
           if match:
              mod_name = match.group(1)
              match = mod_re.match(mod_name)
              if match:
                 mod = match.expand("\\1.mod")
                 if mod not in mods:
                    mods.append(mod)
    f.close()
#
# Returns the list of dependencies for each file
#
    return (deps, mods)

#
# Writes the output file
#

def write_deps(outf, filename, deps, mods):
    filebase, fileext = os.path.splitext(filename)
    outf.write("%s.o: %s\n" % (filebase, " ".join(deps)))
    for mod in mods:
        outf.write("%s: %s\n" % (mod, " ".join(deps)))

#
# Main
#
def main():

#
# Gets the files
#
    opts,filenames = getopt.getopt(sys.argv[1:],"",[])
#
# Outputs on Standard Output
#
    outf = sys.stdout
#
# Header
#
    outf.write("# DO NOT EDIT --- auto-generated file\n")
#
# Loops over the files
#
    for filename in filenames:
#
# Find the dependencies
#
        (deps, mods) = get_deps_and_mods(filename)
        if deps:
#
# Outputs them
#
           write_deps(outf, filename, deps, mods)
#
# Close Standart Output
#
    outf.close()

if __name__ == "__main__":
    main()
