#!/usr/bin/env python3
import os

# --- 入力 ---
material = input("material: ")
defect = input("defect: ")
method = input("method: ")
target_atom = input("target atom: ")

output_base_dir = os.path.expanduser("~/Library/CloudStorage/Box-Box/output")
output_dir = f"{method}/{material}/{defect}"

full_output_dir = os.path.join(output_base_dir, output_dir)


# ---- export only ----
print(f'export LATTICE_CONST=3.076599647')
print(f'export LATTICE_TYPE="Diamond"')
print(f'export OUTPUT_PATH="{full_output_dir}"')
print(f'export ATOM={target_atom}')
print(f'export Burgers=3.0e-10←←←←←←←←←←←←←←任意に変更')