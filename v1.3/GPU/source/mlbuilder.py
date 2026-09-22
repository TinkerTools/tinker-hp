import argparse
import json
import os
from pathlib import Path
import cffi
import re
import site
import sys
ffibuilder = cffi.FFI()

parser = argparse.ArgumentParser(description="Build the Deep-HP CFFI plugin")
parser.add_argument("--module", help="CFFI Python interface input")
parser.add_argument("--header", help="C interface header input")
parser.add_argument("--output", help="platform-native plugin output path")
parser.add_argument("--python-library-dir", help="directory containing the selected shared libpython")
args = parser.parse_args()
cmake_mode = args.module is not None
if not args.module:
    args.module = "mlinterface.py"
if not args.header:
    args.header = "mlinterface.h"
if not args.output:
    args.output = "libmlinterface.dylib"

with open(args.module) as f:
    module = f.read()

default_model_dir = "." if cmake_mode else f"{os.getcwd()}/../ml_models"
module = re.sub(r'''_default_model_dir *= *["'].*["']''', f'_default_model_dir = "{default_model_dir}"', module)
with open(args.header) as f:
    header = f.read()

ffibuilder.embedding_api(header)
link_args = []
source_prefix = ""
if sys.platform.startswith("linux"):
    # Without a SONAME, ld records the relative build-tree path it receives
    # from CMake (for example lib/libmlinterface.so) in executables' DT_NEEDED
    # entry.  A slash in DT_NEEDED bypasses RUNPATH, making the plugin
    # unresolvable at runtime.  Publish the basename as the ELF SONAME.
    link_args.append("-Wl,-soname,libmlinterface.so")
    if args.python_library_dir:
        link_args.append(f"-Wl,-rpath,{args.python_library_dir}")
    embedded_site_packages = os.pathsep.join(site.getsitepackages())
    source_prefix = f'''
    #include <stdio.h>
    #include <stdlib.h>
    #include <string.h>

    __attribute__((constructor))
    static void tinkerhp_configure_deephp_pythonpath(void)
    {{
        const char *site_packages = {json.dumps(embedded_site_packages)};
        const char *current = getenv("PYTHONPATH");
        if (current == NULL || *current == '\\0') {{
            setenv("PYTHONPATH", site_packages, 1);
            return;
        }}

        size_t length = strlen(site_packages) + 1 + strlen(current) + 1;
        char *pythonpath = malloc(length);
        if (pythonpath == NULL) {{
            return;
        }}
        snprintf(pythonpath, length, "%s:%s", site_packages, current);
        setenv("PYTHONPATH", pythonpath, 1);
        free(pythonpath);
    }}
    '''
elif sys.platform == "darwin" and args.python_library_dir:
    link_args.append(f"-Wl,-rpath,{args.python_library_dir}")
ffibuilder.set_source("mlplugin", f'''
    {source_prefix}
    #include "{Path(args.header).resolve()}"
''', extra_link_args=link_args)

ffibuilder.embedding_init_code(module)
ffibuilder.compile(target=str(Path(args.output).resolve()), verbose=True)
