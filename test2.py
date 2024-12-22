import time
from delaunay import sample_integer_points, Triangulation
import matplotlib.pyplot as plt

# Parameters
sizes = range(100, 2100, 100)  # Number of points to test
runtime_with_tree = []   # For use_tree=True
runtime_without_tree = []  # For use_tree=False

# Measure runtime for each n
for n in sizes:
    # Generate points
    points = sample_integer_points(n)
    
    # Case 1: Using SegmentTree (use_tree=True)
    start_time = time.time()
    T_with_tree = Triangulation(points, use_tree=True, make_legal=False)
    T_with_tree.random_incremental()
    runtime_with_tree.append(time.time() - start_time)
    
    # Case 2: Without SegmentTree (use_tree=False)
    start_time = time.time()
    T_without_tree = Triangulation(points, use_tree=False, make_legal=False)
    T_without_tree.random_incremental()
    runtime_without_tree.append(time.time() - start_time)

print("with tree", runtime_with_tree, "\n")
print("no tree", runtime_without_tree)

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(sizes, runtime_with_tree, label='With SegmentTree', marker='o')
plt.plot(sizes, runtime_without_tree, label='No SegmentTree', marker='x')
plt.xlabel('n')
plt.ylabel('Runtime / s')
plt.legend()
plt.grid(True)
plt.show()
