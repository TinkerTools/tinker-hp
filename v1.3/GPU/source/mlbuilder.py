import cffi
import os
import re
ffibuilder = cffi.FFI()

modulename="mlinterface"

with open(modulename+'.py','r') as f:
    module = f.read()

module=re.sub('''_default_model_dir *= *["'].*["']''',f'_default_model_dir = "{os.getcwd()}/../ml_models"',module)
with open(modulename+".h", "r") as f:
    header = f.read()

ffibuilder.embedding_api(header)
ffibuilder.set_source("mlplugin", f'''
    #include "{modulename}.h"
''')

ffibuilder.embedding_init_code(module)
ffibuilder.compile(target=f"lib{modulename}.dylib", verbose=True)
