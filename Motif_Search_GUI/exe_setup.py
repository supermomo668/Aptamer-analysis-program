from cx_Freeze import setup, Executable
import os

base = None

executables = [Executable("AptamerApp.py", base=base)]

packages = ["idna"]
options = {
    'build_exe': {
        'packages':packages
    },
}
## Add Environment
os.environ['TCL_LIBRARY'] = r'C:\ProgramData\Anaconda3\tcl\tcl8.6'
os.environ['TK_LIBRARY'] = r'C:\ProgramData\Anaconda3\tcl\tk8.6'

setup(
    name = "Aptamer Integrate",
    options = options,
    version = "0.9",
    description = 'Aptamer Integral Applications',
    executables = executables
)