import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def parse_path_file(filepath):
    lines = Path(filepath).read_text().strip().split('\n')
    
    num_vectors = None
    vector_size = None
    data_lines  = []
    
    for line in lines:
        if line.startswith('# num_vectors'):
            num_vectors = int(line.split()[-1])
        elif line.startswith('# vector_size'):
            vector_size = int(line.split()[-1])
        elif line.startswith('#'):
            continue
        else:
            data_lines.append(line)
    
    # Parse all numbers
    all_numbers = []
    for line in data_lines:
        all_numbers.extend([float(x) for x in line.split()])
    
    # Split into configurations
    configs = []
    for i in range(num_vectors):
        start = i * vector_size
        vec   = np.array(all_numbers[start:start + vector_size])
        
        n_vertices = (vector_size - 1) // 4  # 3n + n + 1 = 4n + 1
        
        positions = vec[:3 * n_vertices].reshape(n_vertices, 3)
        thetas    = vec[3 * n_vertices:3 * n_vertices + n_vertices]
        Theta     = vec[-1]
        
        configs.append({
            "positions": positions,
            "thetas":    thetas,
            "Theta":     Theta,
        })
    
    return configs

def analyze_twist(configs):
    n = len(configs)
    
    print(f"\n{'='*60}")
    print(f"{'Config':<8} {'Theta':>10} {'mean(theta)':>12} "
          f"{'std(theta)':>12} {'max|theta|':>12}")
    print(f"{'='*60}")
    
    Thetas     = []
    mean_theta = []
    std_theta  = []
    max_theta  = []
    
    for i, cfg in enumerate(configs):
        T  = cfg["Theta"]
        th = cfg["thetas"]
        Thetas.append(T)
        mean_theta.append(np.mean(th))
        std_theta.append(np.std(th))
        max_theta.append(np.max(np.abs(th)))
        print(f"{i:<8} {T:>10.4f} {np.mean(th):>12.4f} "
              f"{np.std(th):>12.4f} {np.max(np.abs(th)):>12.4f}")
    
    print(f"{'='*60}\n")
    
    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    
    x = np.arange(n)
    
    axes[0,0].plot(x, Thetas, marker='o', color='#1f5e8c', linewidth=2)
    axes[0,0].set_title('Global twist $\\Theta$ along path')
    axes[0,0].set_xlabel('Configuration')
    axes[0,0].set_ylabel('$\\Theta$')
    axes[0,0].grid(True, alpha=0.3)
    
    axes[0,1].plot(x, mean_theta, marker='o', color='#c0392b', linewidth=2)
    axes[0,1].set_title('Mean $\\theta_i$ along path')
    axes[0,1].set_xlabel('Configuration')
    axes[0,1].set_ylabel('mean$(\\theta_i)$')
    axes[0,1].grid(True, alpha=0.3)
    
    axes[1,0].plot(x, std_theta, marker='o', color='#27964e', linewidth=2)
    axes[1,0].set_title('Std of $\\theta_i$ along path')
    axes[1,0].set_xlabel('Configuration')
    axes[1,0].set_ylabel('std$(\\theta_i)$')
    axes[1,0].grid(True, alpha=0.3)
    
    # theta profile for each config
    for i, cfg in enumerate(configs):
        axes[1,1].plot(cfg["thetas"], 
                       label=f'K_{i}', 
                       alpha=0.7, linewidth=1.5)
    axes[1,1].set_title('$\\theta_i$ profile per configuration')
    axes[1,1].set_xlabel('Vertex $i$')
    axes[1,1].set_ylabel('$\\theta_i$')
    axes[1,1].legend(fontsize=7)
    axes[1,1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('twist_analysis.png', dpi=150, bbox_inches='tight')
    print("Saved: twist_analysis.png")
    plt.show()

if __name__ == "__main__":
    import sys
    filepath = sys.argv[1] if len(sys.argv) > 1 else "path.txt"
    configs  = parse_path_file(filepath)
    print(f"Loaded {len(configs)} configurations, "
          f"{configs[0]['positions'].shape[0]} vertices each.")
    analyze_twist(configs)
