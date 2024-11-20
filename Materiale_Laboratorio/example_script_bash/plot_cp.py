import matplotlib.pyplot as plt

# Function to read two columns of points from a file
def read_points(filename):
    x = []
    y = []
    
    with open(filename, 'r') as file:
        next(file)  # Skip the first line
        for line in file:
            # Split the line into two values
            values = line.split()
            # Append the values to x and y lists
            x.append(float(values[0]))
            y.append(-float(values[1]))
    
    return x, y

# Example usage
x1, y1 = read_points('naca_0006_1_cp.txt')
x2, y2 = read_points('naca_0006_2_cp.txt')

plt.scatter(x1, y1, color='blue', marker='o')    
plt.scatter(x2, y2, color='red', marker='o')
plt.title('Test cp')
plt.xlabel('X values')
plt.ylabel('-Cp values')
plt.grid(True)
plt.show()
