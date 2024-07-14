import math

class Circle:
    """Circle class represents a circle with the center of x,y"""
    
    def __init__(self, x, y, r):
        """Create the center and radius of the circle"""
        self.x = x
        self.y = y
        self.radius = r
    
    def perimeter(self):
        """Compute the perimeter of the circle"""
        return 2 * math.pi * self.radius
    
    def area(self):
        """Compute the area of the circle"""
        return math.pi * (self.radius ** 2)

class Point:
    """Point class represents a point in 2D space"""
    
    def __init__(self, x, y):
        """Create the x and y coordinates of the point"""
        self.x = x
        self.y = y

    def shortest_distance_to_circle(self, circle):
        """Calculate the shortest distance from the Point to the circumference of a Circle"""
        if not isinstance(circle, Circle):
            raise ValueError("The argument must be an instance of Circle")
        
        # Distance from the point to the circle's center
        distance_to_center = math.sqrt((self.x - circle.x) ** 2 + (self.y - circle.y) ** 2)
        
        # The shortest distance to the circumference of the circle is the absolute value of 
        # the distance from the point to the circle's center minus the circle's radius
        distance_to_circumference = abs(distance_to_center - circle.radius)
        
        # Return the shortest distance rounded to 2 decimal places
        return round(distance_to_circumference, 2)

# Test the code with Circle(-3, 3, 5) and Point(-2, 0)
circle = Circle(-3, 3, 5)
point = Point(-2, 0)
print(point.shortest_distance_to_circle(circle))  # Should print the shortest distance to the circle
