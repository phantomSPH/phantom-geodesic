"""
Python function that executes bash commands.
"""
import subprocess

def BASH(command):
   return subprocess.check_output(command,shell=True).decode().strip()
