import subprocess
from os import listdir, system

system("gcc -o intopt intopt.c")

for _dir in listdir("test_cases"):
    command = f"./intopt < test_cases/{_dir}/i"
    correct = open(f"test_cases/{_dir}/sol", "r").read()

    output = subprocess.check_output(command, shell=True, text=True).strip()
    print(f"Input: {_dir}")
    print(f"Is correct: {output == correct}")
    print()
