import sys
import os
import runpy

script_dir = os.path.dirname(os.path.abspath(__file__))

if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

#print(f"Running package 'genecaller' found in: {script_dir}")
runpy.run_module('genecaller', run_name='__main__', alter_sys=True)

