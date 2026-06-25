import matplotlib.pyplot as plt
import numpy as np
import re
import sys
from pathlib import Path

def parse_neb_file(filepath):
    text = Path(filepath).read_text()
    
    # Max NEB Force als Key
    max_force_match = re.search(r'Max NEB Force:\s*([\d.e+\-]+)', text)
    iterations_match = re.search(r'after\s+([\d]+)\s+iterations', text)
    time_match = re.search(r'took:\s*([\d]+)\s*ms', text)
    
    def parse_list(label):
        m = re.search(label + r'[\s\S]*?\n([\d., e+\-]+)', text)
        if m:
            return [float(x) for x in re.findall(r'[\d.e+\-]+', m.group(1))]
        return []

    return {
        "max_force":  float(max_force_match.group(1)) if max_force_match else None,
        "iterations": int(iterations_match.group(1)) if iterations_match else None,
        "time_ms":    int(time_match.group(1)) if time_match else None,
        "energies":   parse_list("Energies"),
        "neb_forces": parse_list("NEB Forces"),
        "gradients":  parse_list("Gradients"),
    }

def plot_snapshots(filepaths):
    snapshots = {}
    for fp in filepaths:
        data = parse_neb_file(fp)
        if data["max_force"] is not None:
            snapshots[data["max_force"]] = data
            print(f"Loaded {fp}: max_force={data['max_force']:.5f}, "
                  f"it={data['iterations']}, time={data['time_ms']}ms")

    if not snapshots:
        print("Keine Daten gefunden.")
        return

    n_snaps = len(snapshots)
    colors  = plt.cm.viridis(np.linspace(0, 1, max(n_snaps, 2)))

    fig, axes = plt.subplots(3, 1, figsize=(9, 11))
    titles = ["Energy",    "NEB Force",    "Gradient"]
    keys   = ["energies",  "neb_forces",   "gradients"]

    for ax, title, key in zip(axes, titles, keys):
        for idx, (max_f, data) in enumerate(sorted(snapshots.items(), reverse=True)):
            x = np.arange(len(data[key]))
            label = f"max_F={max_f:.5f}"
            if data["iterations"]:
                label += f"  it={data['iterations']}"
            ax.plot(x, data[key], marker='o', label=label,
                    color=colors[idx], linewidth=2, markersize=5)
        ax.set_title(title)
        ax.set_xlabel("Image")
        ax.set_ylabel(title)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("neb_progress.png", dpi=150)
    print("Gespeichert: neb_progress.png")
    plt.show()

# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    if len(sys.argv) > 1:
        files = sys.argv[1:]
    else:
        # alle .txt Dateien im aktuellen Verzeichnis
        files = sorted(Path(".").glob("*.txt"))

    plot_snapshots([str(f) for f in files])
