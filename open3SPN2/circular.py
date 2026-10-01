"""Geometry and topology helpers for circular DNA (plasmids, minicircles).

A circular chain is a chain whose last residue is bonded to its first one. Open3SPN2 keeps the
residues of a circular chain numbered consecutively (``start`` ... ``start + length - 1``) and wraps
every residue offset of the bonded terms around the chain, so the residue after the last one is the
first one again. The functions here provide that wrapping and the tools to build a circular
structure from a straight duplex:

* :func:`circular_chain_info` and :func:`residue_wrapper` describe and wrap the circular chains.
* :class:`CircularBend` maps a duplex onto a planar circle, adding or removing turns of the helix so
  the two ends join with a chosen linking number (relaxed, under- or overwound).
* :func:`linking_number` measures the linking number of two closed curves, e.g. the two strands of a
  plasmid along a trajectory.
"""
import numpy as np


def circular_chain_info(atoms, circular_chains):
    """Returns ``{chainID: (first resSeq, number of residues)}`` for every circular chain.

    The residues of a circular chain must be numbered consecutively, as the residue that follows
    the last one is computed by wrapping the residue number around the chain."""
    info = {}
    for chain in circular_chains:
        residues = np.unique(atoms.loc[atoms['chainID'] == chain, 'resSeq'].astype(int))
        if len(residues) == 0:
            raise ValueError(f'Circular chain {chain} is not in the structure')
        start, length = int(residues.min()), int(residues.max() - residues.min() + 1)
        if len(residues) != length:
            raise ValueError(f'The residues of the circular chain {chain} must be numbered consecutively')
        info[chain] = (start, length)
    return info


def residue_wrapper(atoms, circular_chains):
    """Returns a function ``wrap(chainID, resSeq)`` that maps a residue number past the end (or
    before the start) of a circular chain back into the chain. Linear chains are not changed."""
    info = circular_chain_info(atoms, circular_chains)

    def wrap(chain, resSeq):
        if chain in info:
            start, length = info[chain]
            return start + (resSeq - start) % length
        return resSeq

    return wrap


def kabsch(mobile, reference):
    """Optimal rotation ``R`` and translation ``t`` so that ``mobile @ R.T + t`` fits `reference`."""
    mobile = np.asarray(mobile, dtype=float)
    reference = np.asarray(reference, dtype=float)
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    h = (mobile - mobile_center).T @ (reference - reference_center)
    u, s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ np.diag([1, 1, d]) @ u.T
    translation = reference_center - mobile_center @ rotation.T
    return rotation, translation


class CircularBend:
    """Maps a DNA duplex onto a planar circle, one rigid base pair at a time.

    Each base pair gets a local frame: its origin and tangent follow the helical axis (a line fitted
    to the base pair centers over about one helical turn) and its first axis points from one sugar to
    the other. The base pairs are then placed rigidly along a circle lying on the xy plane, centered at the
    origin, keeping the rise and the twist between consecutive base pairs. Any intrinsic curvature of
    the duplex is replaced by the uniform curvature of the circle. An extra twist, distributed
    uniformly along the circle, makes the end of the duplex join its start with an integer linking
    number.

    Parameters
    ----------
    centers : (n, 3) array
        Center of each base pair, ordered along the duplex.
    references : (n, 3) array
        A vector fixed to each base pair, from the sugar of the first strand to the sugar of the
        second one, used to measure the twist.
    closed : bool
        If True, base pair `last` is a copy of base pair `first` translated one step along the helix,
        so the length and twist of the step that closes the circle are known exactly. If False, the
        closing step is the average step of the duplex.
    delta_linking_number : int
        Turns added to (positive, overwound) or removed from (negative, underwound) the relaxed
        linking number Lk0, the number of helical turns of the duplex rounded to an integer.
    superhelical_density : float, optional
        Alternative to `delta_linking_number`; the change is ``round(sigma * Lk0)``.
    first, last : int, optional
        Base pairs that form the circle (by default all of them). Base pairs outside this range are
        still mapped, continuing the circle, which is useful to pad the ends of the duplex so all the
        base pairs of the circle have a complete neighborhood to compute their frames.
    window : int
        Number of base pairs used to fit the helical axis.
    """

    def __init__(self, centers, references, closed=False, delta_linking_number=0, superhelical_density=None,
                 first=0, last=None, window=11):
        centers = np.asarray(centers, dtype=float)
        references = np.asarray(references, dtype=float)
        n = len(centers)
        last = n - 1 if last is None else last
        if last - first < 2:
            raise ValueError('At least three base pairs are needed to bend a duplex into a circle')

        # Local frames of the base pairs: the helical axis is the line that best fits the centers of the
        # `window` base pairs around each base pair (shifted inwards at the ends of the duplex)
        window = min(window, n)
        origins = np.zeros((n, 3))
        tangents = np.zeros((n, 3))
        for i in range(n):
            start = min(max(i - window // 2, 0), n - window)
            points = centers[start:start + window]
            mean = points.mean(axis=0)
            direction = np.linalg.svd(points - mean)[2][0]
            if (points[-1] - points[0]) @ direction < 0:
                direction = -direction
            origins[i] = mean + ((centers[i] - mean) @ direction) * direction
            tangents[i] = direction
        e1 = references - np.sum(references * tangents, axis=1, keepdims=True) * tangents
        e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
        e2 = np.cross(tangents, e1)
        self.origins = origins
        self.frames = np.stack([e1, e2, tangents], axis=2)  # columns e1, e2, t

        # Rise and twist of each step; the twist is the rotation of e1 around the axis after moving it
        # along the axis without twisting it (parallel transport)
        rise = np.linalg.norm(np.diff(origins, axis=0), axis=1)
        twist = np.zeros(n - 1)
        for i in range(n - 1):
            moved = _minimal_rotation(tangents[i], tangents[i + 1]) @ e1[i]
            twist[i] = np.arctan2(np.cross(moved, e1[i + 1]) @ tangents[i + 1], moved @ e1[i + 1])
        contour = np.concatenate([[0], np.cumsum(rise)])
        turns = np.concatenate([[0], np.cumsum(twist)])
        self.contour = contour - contour[first]
        self.turns = turns - turns[first]

        if closed:
            self.circumference = self.contour[last]
            self.twist = self.turns[last]  # radians
        else:
            self.circumference = self.contour[last] * (last - first + 1) / (last - first)
            self.twist = self.turns[last] * (last - first + 1) / (last - first)
        self.radius = self.circumference / (2 * np.pi)

        self.relaxed_linking_number = int(np.round(self.twist / (2 * np.pi)))
        if superhelical_density is not None:
            delta_linking_number = int(np.round(superhelical_density * self.relaxed_linking_number))
        self.delta_linking_number = int(delta_linking_number)
        self.linking_number = self.relaxed_linking_number + self.delta_linking_number
        self.extra_twist = 2 * np.pi * self.linking_number - self.twist

        # Target frames on the circle; e1 -> z, e2 -> radial direction, t -> tangent of the circle when untwisted
        theta = 2 * np.pi * self.contour / self.circumference
        alpha = self.turns + self.extra_twist * self.contour / self.circumference
        z = np.array([0, 0, 1.])
        radial = np.stack([np.cos(theta), np.sin(theta), np.zeros(n)], axis=1)
        tangent = np.stack([-np.sin(theta), np.cos(theta), np.zeros(n)], axis=1)
        target_e1 = np.cos(alpha)[:, None] * z + np.sin(alpha)[:, None] * radial
        target_e2 = -np.sin(alpha)[:, None] * z + np.cos(alpha)[:, None] * radial
        self.target_origins = self.radius * radial
        self.target_frames = np.stack([target_e1, target_e2, tangent], axis=2)

    def apply(self, xyz, base_pair):
        """Returns the bent coordinates of `xyz` (same units as the base pair centers). `base_pair` is the
        index of the base pair each coordinate belongs to."""
        xyz = np.asarray(xyz, dtype=float)
        base_pair = np.asarray(base_pair, dtype=int)
        local = np.einsum('nji,nj->ni', self.frames[base_pair], xyz - self.origins[base_pair])
        return np.einsum('nij,nj->ni', self.target_frames[base_pair], local) + self.target_origins[base_pair]


def _minimal_rotation(a, b):
    """Rotation matrix that takes the unit vector `a` to the unit vector `b` around their common normal"""
    v = np.cross(a, b)
    c = a @ b
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx / (1 + c)


def linking_number(curve1, curve2):
    """Gauss linking number of two closed polygonal curves (``(n, 3)`` arrays, the last point is
    connected back to the first one).

    For the two strands of a circular duplex, e.g. the sugar beads of each strand, this is the
    linking number Lk = Tw + Wr of the plasmid, which is constant along an unbroken simulation. It is
    positive for a right handed double helix. Uses the exact solid angle of every pair of segments
    (Klenin & Langowski, Biopolymers 54, 307 (2000))."""
    p1 = np.asarray(curve1, dtype=float)
    p3 = np.asarray(curve2, dtype=float)
    p2 = np.roll(p1, -1, axis=0)
    p4 = np.roll(p3, -1, axis=0)
    p1, p2 = p1[:, None], p2[:, None]
    p3, p4 = p3[None, :], p4[None, :]

    r13, r14, r23, r24 = p3 - p1, p4 - p1, p3 - p2, p4 - p2

    def normal(a, b):
        n = np.cross(a, b)
        norm = np.linalg.norm(n, axis=-1, keepdims=True)
        return np.divide(n, norm, out=np.zeros_like(n), where=norm > 0)

    n1, n2, n3, n4 = normal(r13, r14), normal(r14, r24), normal(r24, r23), normal(r23, r13)

    def asin(x):
        return np.arcsin(np.clip(x, -1, 1))

    omega = (asin(np.sum(n1 * n2, -1)) + asin(np.sum(n2 * n3, -1)) +
             asin(np.sum(n3 * n4, -1)) + asin(np.sum(n4 * n1, -1)))
    sign = np.sign(np.sum(np.cross(p4 - p3, p2 - p1) * r13, -1))
    return float(-np.sum(omega * sign) / (4 * np.pi))
