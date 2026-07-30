import os, re

folder = os.path.dirname(os.path.abspath(__file__))

N = 0
for fname in sorted(f for f in os.listdir(folder) if re.match(r'\w+\.txt', f)):
    with open(os.path.join(folder, fname)) as f:
        n = sum(1 for line in f if line.strip() and line[0].isdigit())
        N = N + n
    print(f"{fname}: {n} data points")
print(N)