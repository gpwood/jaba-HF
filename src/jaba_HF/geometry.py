import numpy as np
import math
from dataclasses import dataclass, field
from typing import List

@dataclass
class point_xyz:
    """ class to store an xyz point with atom """
    atom: str
    x: float
    y: float
    z: float

@dataclass
class geometry:
    """ a collection of point xyz's """
    points: List[point_xyz] = field(default_factory=list)
     
    def np_matrix():
        matrix_np = []
        for point in points:
            matrix_np.append(np.array(point.x,point.y,point.z))
        return np.array(matrix_np)

    def append(point):
        points.append(point)
