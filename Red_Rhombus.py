# Adjusting the code to make spikes face outward and perpendicular to the legs
import numpy as np
import matplotlib.pyplot as plt
plt.savefig('C:\\Users\\jaan\\OneDrive\\Pictures\\Saved Pictures\\red_rhombus.png', bbox_inches='tight', pad_inches=0, dpi=100, facecolor='black')
ax = plt.gca()
ax.set_facecolor('black') # Ensure the background is black

# Re-defining the coordinates of the rhombus
x = np.array([0, 1, 0, -1, 0])
y = np.array([2, 0, -2, 0, 2])

# Calculating slopes of legs to determine perpendicular direction for spikes
slope_legs = np.diff(y) / np.diff(x)

# Handling infinite slope (vertical lines)
slope_legs = np.append(slope_legs, -1/slope_legs[1])  # Use the reciprocal of the slope for perpendicular direction

# Plotting the rhombus
plt.plot(x, y, 'r-', linewidth=5) # 'r-' for red line

# Adding spikes: they should be perpendicular, so we'll adjust directions accordingly
for i in range(4):
    mx, my = (x[i] + x[i+1]) / 2, (y[i] + y[i+1]) / 2  # Midpoints
    if np.isinf(slope_legs[i]):  # Vertical leg
        spike_direction = np.array([[0, 1/3], [0, -1/3]])
    else:
        dx = 1/3 / np.sqrt(1 + slope_legs[i]**2)
        dy = slope_legs[i] * dx
        spike_direction = np.array([[-dy, dx], [dy, -dx]])

    # Plotting two spikes per side to ensure they are outward and perpendicular
    for dir_vector in spike_direction:
        spike_end_x = mx + dir_vector[0]
        spike_end_y = my + dir_vector[1]
        plt.plot([mx, spike_end_x], [my, spike_end_y], 'r-', linewidth=5)

plt.axis('equal')
plt.axis('off')
plt.gca().set_facecolor('black')
plt.savefig('C:/Users/jaan/OneDrive/Pictures/Saved Pictures/red_rhombus.png', facecolor=plt.gcf().get_facecolor(), bbox_inches='tight', pad_inches=0, dpi=100)
plt.show()
