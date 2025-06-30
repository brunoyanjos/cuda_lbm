import pandas as pd
import numpy as py
import matplotlib.pyplot as plt

tke_path = 'LDC/001/data10002500.f'

with open(tke_path, 'rb') as file:  # 'rb' = read binary
    raw_data = file.read()
    print(f"Read {len(raw_data)} bytes")

