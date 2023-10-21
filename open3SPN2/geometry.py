import numpy as np

class Transform:
    def __init__(self, dx=0, dy=0, dz=0, rx=0, ry=0, rz=0):
        """
        Initialize the Transform object with displacements and rotations.
        
        Parameters:
        - dx, dy, dz: Displacements in the x, y, z directions respectively (Angstroms).
        - rx, ry, rz: Rotations around the x, y, z axes respectively (Radians).
        """
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.rx = rx
        self.ry = ry
        self.rz = rz
        
        # Angle needed to rotate around the z axis to achieve a single rotation around the y axis
        self.phi = np.arctan2(self.rx, self.ry)
        
        # Magnitude of the combined x and y rotation
        self.yx_ang = np.sqrt(self.rx**2 + self.ry**2)

    @property
    def full_rotation(self):
        """
        Compute the full rotation transformation matrix that transitions from one state to another.
        
        Returns:
        - A 3x3 rotation matrix.
        """
        rotz1 = self.rotz(0.5 * self.rz - self.phi)
        rotyx = self.roty(self.yx_ang)
        rotz2 = self.rotz(0.5 * self.rz + self.phi)
        return rotz1 @ rotyx @ rotz2

    @property
    def half_rotation(self):
        """
        Compute the rotation matrix for half the rotation.
        
        This aims to rotate two objects equally relative to a reference frame.
        
        Returns:
        - A 3x3 rotation matrix.
        """
        rotz1 = self.rotz(0.5 * self.rz - self.phi)
        rotyx = self.roty(0.5 * self.yx_ang)
        rotz2 = self.rotz(self.phi)
        return rotz1 @ rotyx @ rotz2

    @property
    def full_displacement(self):
        """
        Compute the displacement vector for transitioning from one state to another.
        
        Returns:
        - A 1x3 displacement vector.
        """
        displacement = np.array([self.dx, self.dy, self.dz])
        return np.dot(displacement, self.half_rotation.T)

    @property
    def half_displacement(self):
        """
        Compute the midpoint of the displacement vector.
        
        Returns:
        - A 1x3 displacement vector.
        """
        return 0.5 * self.full_displacement
    
    @staticmethod
    def rotz(theta):
        """
        Compute the rotation matrix around the Z-axis.
        
        Parameters:
        - theta: Rotation angle in radians.
        
        Returns:
        - A 3x3 rotation matrix.
        """
        c, s = np.cos(theta), np.sin(theta)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

    @staticmethod
    def roty(theta):
        """
        Compute the rotation matrix around the Y-axis.
        
        Parameters:
        - theta: Rotation angle in radians.
        
        Returns:
        - A 3x3 rotation matrix.
        """
        c, s = np.cos(theta), np.sin(theta)
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def test_default_values():
    transform = Transform()
    
    np.testing.assert_almost_equal(transform.full_rotation, np.eye(3), decimal=7)
    np.testing.assert_almost_equal(transform.half_rotation, np.eye(3), decimal=7)
    np.testing.assert_almost_equal(transform.full_displacement, [0, 0, 0], decimal=7)
    np.testing.assert_almost_equal(transform.half_displacement, [0, 0, 0], decimal=7)

def test_given_values():
    transform = Transform(0.07, -0.19, 0.07, 
                          1.8/180*np.pi, -15.0/180*np.pi, 1.5/180*np.pi)
    
    expected_full_rotation = np.array([[0.965592, -0.029813, -0.258348], 
                                       [0.021636,  0.999173, -0.034438], 
                                       [0.259161,  0.027663,  0.965438]])
    expected_half_rotation = np.array([[0.991374, -0.014114, -0.130305],
                                       [0.011951,  0.999778, -0.017370],
                                       [0.130521,  0.015662,  0.991322]])
    expected_full_displacement = [0.062957, -0.190337, 0.075553]
    expected_half_displacement = [0.031478, -0.095169, 0.037777]
    
    np.testing.assert_almost_equal(transform.full_rotation, expected_full_rotation, decimal=6)
    np.testing.assert_almost_equal(transform.half_rotation, expected_half_rotation, decimal=6)
    np.testing.assert_almost_equal(transform.full_displacement, expected_full_displacement, decimal=6)
    np.testing.assert_almost_equal(transform.half_displacement, expected_half_displacement, decimal=6)

if __name__=='__main__':
    test_default_values()
    test_given_values()