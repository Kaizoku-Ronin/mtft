"""The inputs of the comparison between ``mtft.plane`` and the Modular Plane lab's own core (version 8.3): matrices for the kernel,
sites of the quantity lattice for the prospector.  Deterministic (seeded); writes inputs.json beside this file.

    python make_inputs.py
"""
import json
import random
from pathlib import Path

NAMED = {
    "pell": [[12, 1], [143, 12]], "cat": [[2, 1], [1, 1]], "dual unit": [[0, 1], [0, 0]], "S": [[0, -1], [1, 0]], "T": [[1, 1], [0, 1]],
    "s10": [[84, -47], [143, -80]], "order 3": [[0, 1], [-1, -1]], "order 6": [[1, -1], [1, 0]], "identity": [[1, 0], [0, 1]],
    "minus identity": [[-1, 0], [0, -1]], "golden, det -1": [[1, 1], [1, 0]], "U V^3": [[4, 1], [3, 1]],
}
SITES = [[-1, -1, 4, 2], [1, -1, -1, 0], [0, 0, 0, 0], [1, 2, -2, 0], [0, -1, 0, 0], [-1, -2, 4, 2], [1, -1, -3, 0], [0, -1, 2, 0], [1, 1, -3, -1]]


def inputs(n_matrices: int = 40, n_sites: int = 40) -> dict:
    rng = random.Random(83)
    matrices = [m for m in NAMED.values()]
    for _ in range(n_matrices):
        scale = 10 ** rng.uniform(-1, 2)
        matrices.append([[rng.uniform(-1, 1) * scale, rng.uniform(-1, 1) * scale], [rng.uniform(-1, 1) * scale, rng.uniform(-1, 1) * scale]])
    rng = random.Random(7)
    sites = [list(v) for v in SITES] + [[rng.randint(-2, 2), rng.randint(-3, 3), rng.randint(-4, 4), rng.randint(-2, 2)] for _ in range(n_sites)]
    return {"names": list(NAMED), "matrices": matrices, "sites": sites, "t": 0.7}


if __name__ == "__main__":
    out = Path(__file__).resolve().parent / "inputs.json"
    out.write_text(json.dumps(inputs(), indent=None) + "\n", encoding="utf-8")
    print("wrote", out.name)
