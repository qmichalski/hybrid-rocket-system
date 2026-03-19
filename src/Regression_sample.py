import math
import matplotlib.pyplot as plt

def intersection_point(x1, y1, x2, y2, x3, y3, x4, y4):
    # Calculate the coefficients of the line equations
    a1 = y2 - y1
    b1 = x1 - x2
    c1 = a1 * x1 + b1 * y1

    a2 = y4 - y3
    b2 = x3 - x4
    c2 = a2 * x3 + b2 * y3

    # Calculate the determinant
    det = a1 * b2 - a2 * b1

    if det == 0:
        return None  # Lines are parallel

    # Calculate the intersection point
    x = (b2 * c1 - b1 * c2) / det
    y = (a1 * c2 - a2 * c1) / det

    return (x, y)

def point_lies_on_segment(x, y, x1, y1, x2, y2):
    return (min(x1, x2) <= x <= max(x1, x2)) and (min(y1, y2) <= y <= max(y1, y2))

def is_intersection_on_segments(x1, y1, x2, y2, x3, y3, x4, y4):
    intersection = intersection_point(x1, y1, x2, y2, x3, y3, x4, y4)
    if intersection is None:
        return False  # Lines are parallel

    x, y = intersection
    return point_lies_on_segment(x, y, x1, y1, x2, y2) and point_lies_on_segment(x, y, x3, y3, x4, y4)


# finding a parallel line segment shifted by a distance d away from the original line segment away from the origin
def parallel_line_segment(x1, y1, x2, y2, d):
    # Calculate the direction vector of the line segment
    dx = x2 - x1
    dy = y2 - y1
    length = math.sqrt(dx**2 + dy**2)

    if length == 0:
        raise ValueError("The endpoints cannot be the same")

    # Normalize the direction vector
    unit_dx = dx / length
    unit_dy = dy / length

    # Calculate the perpendicular vector (rotated 90 degrees)
    perp_dx = -unit_dy
    perp_dy = unit_dx

    # Scale the perpendicular vector by distance d
    shift_x = perp_dx * d
    shift_y = perp_dy * d

    # Calculate the new endpoints of the parallel line segment
    new_x1 = x1 + shift_x
    new_y1 = y1 + shift_y
    new_x2 = x2 + shift_x
    new_y2 = y2 + shift_y

    return (new_x1, new_y1), (new_x2, new_y2)

def three_point_panel(x1, y1, x2, y2, x3, y3, d=1):
    X1 = x1
    Y1 = y1
    X2 = x2
    Y2 = y2
    X3 = x2
    Y3 = y2
    X4 = x3
    Y4 = y3
    (new_x1, new_y1), (new_x2, new_y2) = parallel_line_segment(X1, Y1, X2, Y2, d)
    (new_x3, new_y3), (new_x4, new_y4) = parallel_line_segment(X3, Y3, X4, Y4, d)

    if is_intersection_on_segments(new_x1, new_y1, new_x2, new_y2, new_x3, new_y3, new_x4, new_y4):
        POI_X, POI_Y = intersection_point(new_x1, new_y1, new_x2, new_y2, new_x3, new_y3, new_x4, new_y4)
        list_of_points = [(new_x1, new_y1), (POI_X, POI_Y), (new_x4, new_y4)]
        decision = True
        print(f"Intersection Point: ({POI_X}, {POI_Y})")
    else:
        list_of_points = [(new_x1, new_y1), (new_x2, new_y2), (new_x3, new_y3), (new_x4, new_y4)]
        decision = False
    return list_of_points, decision

def regression_calculator(geometry_points, regression_step=1):
    post_regression_points = []
    flag = 0
    for i in range(len(geometry_points) - 2):
        decision = False
        x1, y1 = geometry_points[i]
        x2, y2 = geometry_points[i + 1]
        x3, y3 = geometry_points[i + 2]
        panel_points, decision = three_point_panel(x1, y1, x2, y2, x3, y3)
        if i !=0:
            if len(panel_points) == 3:
                if decision:
                    post_regression_points.pop()
                post_regression_points.extend(panel_points[1:3])
            elif len(panel_points) == 4:
                if decision:
                    post_regression_points.pop()
                post_regression_points.extend(panel_points[1:4])
        else:
            post_regression_points.extend(panel_points)
        #post_regression_points.pop()
        
    
    x1, y1 = geometry_points[len(geometry_points)-2]
    x2, y2 = geometry_points[len(geometry_points)-1]
    x3, y3 = geometry_points[0]
    panel_points, decision = three_point_panel(x1, y1, x2, y2, x3, y3)
    if i !=0:
            if len(panel_points) == 3:
                if decision:
                    post_regression_points.pop()
                post_regression_points.extend(panel_points[1:3])
            elif len(panel_points) == 4:
                if decision:
                    post_regression_points.pop()
                post_regression_points.extend(panel_points[1:4])
    else:
        post_regression_points.extend(panel_points)

    # Remove duplicates while preserving order
    seen = set()
    unique_points = []
    for point in post_regression_points:
        if point not in seen:
            seen.add(point)
            unique_points.append(point)
    return unique_points


"""""
original_x, original_y = zip(*geometry_points)
regressed_x, regressed_y = zip(*regressed_points)
plt.scatter(original_x, original_y, color='blue', label='Original Points')
plt.scatter(regressed_x, regressed_y, color='red', label='Regressed Points')
plt.plot(regressed_x + (regressed_x[0],), regressed_y + (regressed_y[0],), color='red', linestyle='-', linewidth=1)
plt.plot(original_x + (original_x[0],), original_y + (original_y[0],), color='blue', linestyle='-', linewidth=1)
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Original vs Regressed Points')
plt.legend()
plt.grid()
plt.show()
"""

def regression_iteration(geometry_points, iterations=7, regression_step=1):
    original_x, original_y = zip(*geometry_points)
    #plt.scatter(original_x, original_y, color='blue', label='Original Points')
    plt.plot(original_x + (original_x[0],), original_y + (original_y[0],), color='blue', linestyle='-', linewidth=1)
    regressed_points = geometry_points
    for i in range(iterations):
        print(f"Iteration {i+1}:")
        regressed_points = regression_calculator(regressed_points, regression_step)
        regressed_x, regressed_y = zip(*regressed_points)
        plt.scatter(regressed_x, regressed_y)#, label=f'Regressed Points - Iteration {i+1}')
        plt.plot(regressed_x + (regressed_x[0],), regressed_y + (regressed_y[0],), linestyle='-', linewidth=1)
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Original vs Regressed Points')
    plt.legend()
    plt.grid()
    plt.show()

a = 20
geometry_points = [(0, a),(4, 4), (a, 0), (4, -4), (0,-a), (-4, -4), (-a, 0), (-4, 4)]
#geometry_points = [(0, a),(4, 4), (a, 0)]

regression_iteration(geometry_points, iterations=15, regression_step=1)
