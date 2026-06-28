import matplotlib.pyplot as plt
import numpy as np
import re
import sys
from pathlib import Path

def parse_neb_file(filepath):
    text = Path(filepath).read_text()
    
    max_force_match    = re.search(r'Max NEB Force:\s*([\d.e+\-]+)', text)
    iterations_match   = re.search(r'after\s+([\d]+)\s+iterations', text)
    time_match         = re.search(r'took:\s*([\d]+)\s*ms', text)
    neb_energie_match  = re.search(r'NEB Energy:\s*([\d.e+\-]+)', text)

    def parse_list(label):
        m = re.search(label + r'[\s\S]*?\n([\d., e+\-]+)', text)
        if m:
            return [float(x) for x in re.findall(r'[\d.e+\-]+', m.group(1))]
        return []
    
    energies        = parse_list("Energies")
    neb_forces      = parse_list("NEB Forces")
    gradients       = parse_list("Gradients")
    spring_forces   = parse_list("Spring Forces")
    bend_energies   = parse_list("Bend Energies")
    twist_energies  = parse_list("Twist Energies")
    stretch_energies = parse_list("Stretch Energies")

    return {
        "max_force":       float(max_force_match.group(1))    if max_force_match    else None,
        "iterations":      int(iterations_match.group(1))     if iterations_match   else None,
        "time_ms":         int(time_match.group(1))           if time_match         else None,
        "E_neb":           float(neb_energie_match.group(1))  if neb_energie_match  else None,
        "energies":        energies,
        "neb_forces":      neb_forces,
        "gradients":       gradients,
        "spring_forces":   spring_forces,
        "bend_energies":   bend_energies,
        "twist_energies":  twist_energies,
        "stretch_energies": stretch_energies,
        # derived summary values
        "max_spring":      max(spring_forces)    if spring_forces    else None,
        "max_E_elastic":   max(energies)         if energies         else None,
        "max_grad":        max(gradients)        if gradients        else None,
    }

def print_table(snapshots):
    """Print a LaTeX-ready summary table to the terminal."""
    header = (
        f"{'max(F_NEB)':<14} "
        f"{'max(F_spring||)':<16} "
        f"{'E_NEB':<14} "
        f"{'max(E_elastic)':<16} "
        f"{'max(grad E)':<14} "
        f"{'Time (ms)':<12} "
        f"{'Iterations'}"
    )
    sep = "-" * len(header)
    print("\n" + sep)
    print(header)
    print(sep)
    for max_f, data in sorted(snapshots.items(), reverse=True):  # largest first
        time_s = f"{data['time_ms']/1000:.0f}" if data['time_ms'] is not None else "N/A"
        print(
            f"{data['max_force']        if data['max_force']    is not None else 'N/A':<14.3g} "
            f"{data['max_spring']       if data['max_spring']   is not None else 'N/A':<16.4g} "
            f"{data['E_neb']            if data['E_neb']        is not None else 'N/A':<14.2g} "
            f"{data['max_E_elastic']    if data['max_E_elastic']is not None else 'N/A':<16.2g} "
            f"{data['max_grad']         if data['max_grad']     is not None else 'N/A':<14.3g} "
            f"{time_s:<12} "
            f"{data['iterations']       if data['iterations']   is not None else 'N/A'}"
        )
    print(sep)

    # LaTeX rows for easy copy-paste
    print("\nLaTeX rows (copy into tabular):")
    print(sep)
    for max_f, data in sorted(snapshots.items(), reverse=True):  # largest first
        def fmt(v): return f"{v:.4g}" if v is not None else "--"
        time_s = f"{data['time_ms']/1000:.0f}" if data['time_ms'] is not None else "--"
        print(
            f"{fmt(data['max_force'])} & "
            f"{fmt(data['max_spring'])} & "
            f"{fmt(data['E_neb'])} & "
            f"{fmt(data['max_E_elastic'])} & "
            f"{fmt(data['max_grad'])} & "
            f"{time_s} s / "
            f"{data['iterations'] if data['iterations'] is not None else '--'} \\\\"
        )
    print(sep + "\n")

def plot_snapshots(filepaths):
    snapshots = {}
    for fp in filepaths:
        data = parse_neb_file(fp)
        if data["max_force"] is not None:
            snapshots[data["max_force"]] = data
            print(f"Loaded {fp}: max_force={data['max_force']:.3g}, "
                  f"it={data['iterations']}, time={data['time_ms']}ms")

    if not snapshots:
        print("Keine Daten gefunden.")
        return

    print_table(snapshots)

    n_snaps = len(snapshots)
    colors  = plt.cm.viridis(np.linspace(0, 1, max(n_snaps, 2)))

    # ---- Plot 1: Main metrics ----
    titles_main = [   "NEB Force",  "Spring Force",  "Energy","Gradient"]
    keys_main   = [ "neb_forces", "spring_forces" ,"energies","gradients"]

    fig1, axes1 = plt.subplots(4, 1, figsize=(9, 11))
    for ax, title, key in zip(axes1, titles_main, keys_main):
        for idx, (max_f, data) in enumerate(sorted(snapshots.items(), reverse=True)):
            vals = data[key]
            if not vals:
                continue
            x     = np.arange(len(vals))
            label = f"max_F={max_f:.3g}"
            if data["iterations"]:
                label += f"  it={data['iterations']}"
            ax.plot(x, vals, marker='o', label=label,
                    color=colors[idx], linewidth=2, markersize=5)
        ax.set_title(title)
        ax.set_xlabel("Configuration")
        ax.set_ylabel(title)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    fig1.tight_layout()
    fig1.savefig("neb_progress.png", dpi=150)
    print("Gespeichert: neb_progress.png")

    # ---- Plot 2: Energy breakdown ----
    titles_energy = ["Bend Energy",   "Twist Energy",   "Stretch Energy"]
    keys_energy   = ["bend_energies", "twist_energies", "stretch_energies"]

    fig2, axes2 = plt.subplots(3, 1, figsize=(9, 8))
    for ax, title, key in zip(axes2, titles_energy, keys_energy):
        for idx, (max_f, data) in enumerate(sorted(snapshots.items(), reverse=True)):
            vals = data[key]
            if not vals:
                continue
            x     = np.arange(len(vals))
            label = f"max_F={max_f:.3g}"
            if data["iterations"]:
                label += f"  it={data['iterations']}"
            ax.plot(x, vals, marker='o', label=label,
                    color=colors[idx], linewidth=2, markersize=5)
        ax.set_title(title)
        ax.set_xlabel("Configuration")
        ax.set_ylabel(title)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    fig2.tight_layout()
    fig2.savefig("energy_breakdown.png", dpi=150)
    print("Gespeichert: energy_breakdown.png")

    plt.show()
# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    if len(sys.argv) > 1:
        files = sys.argv[1:]
    else:
        files = sorted(Path(".").glob("*.txt"))
    plot_snapshots([str(f) for f in files])
