import time
from delaunay import sample_integer_points, Triangulation
import matplotlib.pyplot as plt

# Parameters
sizes = range(100, 2100, 100)  # Number of points to test
runtime_make_legal = []  # For make_legal=True
runtime_postponed = []   # For make_legal=False + naive_delaunay()

# Measure runtime for each n
for n in sizes:
    # Generate points
    points = sample_integer_points(n)
    
    # Case 1: Incremental legalization (make_legal=True)
    start_time = time.time()
    T1 = Triangulation(points, use_tree=False, make_legal=True)
    T1.random_incremental()
    runtime_make_legal.append(time.time() - start_time)
    
    # Case 2: Postponed legalization (make_legal=False + naive_delaunay)
    start_time = time.time()
    T2 = Triangulation(points, use_tree=False, make_legal=False)
    T2.random_incremental()
    T2.naive_delaunay()
    runtime_postponed.append(time.time() - start_time)

print("with increment", runtime_make_legal, "\n")
print("naive", runtime_postponed)

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(sizes, runtime_make_legal, label='Incremental', marker='o')
plt.plot(sizes, runtime_postponed, label='Naive', marker='x')
plt.xlabel('n')
plt.ylabel('Runtime / s')
plt.legend()
plt.grid(True)
plt.show()
