import os
import subprocess
import sys

def BASH(command):
    return subprocess.check_output(command,shell=True).decode().strip()

HOME_DIR = BASH('echo $HOME')
CODE_DIR = HOME_DIR+'/grtest'
CUR_DIR  = BASH('pwd')

try:
    foldername = sys.argv[1]
except:
    quit('No folder name given...')

FOLDER_PATH = CUR_DIR+'/'+foldername

files2copy=[
            CODE_DIR+'/Makefile',
            CUR_DIR+'/splash.*'
            ]

files2link=[
            CODE_DIR+'/scripts/plot.ipynb'
            ]

print(BASH('mkdir -v '+FOLDER_PATH))

for filei in files2copy:
    print(BASH('cp -v '+filei+' '+FOLDER_PATH+'/.'))

for filei in files2link:
    print(BASH('ln -sv '+filei+' '+FOLDER_PATH+'/.'))
