import os
import shutil

homeDir = os.getcwd()

shutil.rmtree(f'{homeDir}/RST_05')
shutil.rmtree(f'{homeDir}/RST_06')
shutil.rmtree(f'{homeDir}/RST_07_snapshots')
shutil.rmtree(f'{homeDir}/TRJ_06')
shutil.rmtree(f'{homeDir}/RST_INIT')
