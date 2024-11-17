import re
import sys

def extract_coordinates(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
    
    # Find the block that starts with `/coor [` and ends with `] def`
    match = re.search(r'/coor\s+\[(.*?)\]\s+def', content, re.DOTALL)
    if not match:
        raise ValueError("Coordinates block not found in the file.")
    
    # Extract the block content
    block_content = match.group(1)
    
    # Find all coordinate pairs using regex
    coordinates = re.findall(r'\[([-\d\.]+)\s+([-\d\.]+)\]', block_content)
    
    # Convert the pairs to a list of tuples (float)
    coordinates = [(float(x), float(y)) for x, y in coordinates]
    return coordinates


def coords2edges(coordinates):
    edges = []
    for i in range(len(coordinates) - 1):
        edges.append((coordinates[i], coordinates[i + 1]))
    return edges

# Check if two edges cross each other
def is_cross(edge1, edge2):
    x1, y1 = edge1[0]
    x2, y2 = edge1[1]
    x3, y3 = edge2[0]
    x4, y4 = edge2[1]
    if (max(x1, x2) < min(x3, x4) or max(x3, x4) < min(x1, x2) or
        max(y1, y2) < min(y3, y4) or max(y3, y4) < min(y1, y2)):
        return False
    if ((x1 - x2) * (y3 - y1) + (y1 - y2) * (x1 - x3)) * ((x1 - x2) * (y4 - y1) + (y1 - y2) * (x1 - x4)) > 0:
        return False
    if ((x3 - x4) * (y1 - y3) + (y3 - y4) * (x3 - x1)) * ((x3 - x4) * (y2 - y3) + (y3 - y4) * (x3 - x2)) > 0:
        return False
    return True

# check if there is cross in a series of edges:
def check_cross(file_path):
    coordinates = extract_coordinates(file_path)
    edges = coords2edges(coordinates)
    for i in range(len(edges)):
        for j in range(i + 2, len(edges)):
            if is_cross(edges[i], edges[j]):
                print(1)
                return True
    print(0)
    return False

if __name__ == "__main__":

    # Usage
    file_path = sys.argv[1]  # Replace with the actual path to your .ps file
    # print(file_path, end=' ')
    check_cross(file_path)
    # coordinates = extract_coordinates(file_path)

    # # Print the coordinates
    # for coord in coordinates:
    #     print(coord)
