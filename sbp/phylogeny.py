"""
Find an orthonormal basis in the Aitchison simplex using a technique
known as sequential binary paritioning.  A phylogenetic tree
will be used as input
"""
from __future__ import division

from skbio import TreeNode
from skbio.stats.composition import clr_inv
import numpy as np

def phylogenetic_basis(treenode):
    """
    Determines the basis based on phylogenetic tree

    Parameters
    ----------
    treenode : skbio.TreeNode
        Phylogenetic tree.  MUST be a bifurcating tree

    """
    basis, _ =  _sequential_binary_partition(treenode)
    return basis


def _sequential_binary_partition(treenode, left_step=0, right_step=0):
    """
    Determines the basis based on phylogenetic tree

    Parameters
    ----------
    treenode : skbio.TreeNode
        Phylogenetic tree.  MUST be a bifurcating tree
    left_step : int
        number of left branches investigated
    right_step : int
        number of right branches investigated

    Returns
    -------
    basis : dict, np.array
        Dictionary of orthonormal bases
    num_leaves: int
        Number of leaves in subtree

    Raises
    ------
    ValueError
        The tree doesn't contain two branches
    ValueError
        The tree doesn't have unique node names

    """
    if len(treenode.children) == 0:
        return {}, 1

    if len(treenode.children) != 2:
        raise ValueError("Not a bifurcating tree!")

    left_basis, r = _sequential_binary_partition(treenode.children[0],
                                                 left_step = left_step + 1,
                                                 right_step = right_step)
    right_basis, s = _sequential_binary_partition(treenode.children[1],
                                                  left_step = left_step,
                                                  right_step = right_step + 1)
    print "Node: %s"%treenode.name,"r=%d,s=%s"%(r, s)
    a = np.sqrt(s / (r*(r+s)))
    b = -1*np.sqrt(r / (s*(r+s)))
    base = clr_inv([0]*right_step + [a]*r + [b]*s + [0]*left_step)
    basis = _merge_two_dicts(left_basis, right_basis)

    basis[treenode.name] = base
    num_leaves = r + s
    return basis, num_leaves

def _merge_two_dicts(x, y):
    '''
    Given two dicts, merge them into a new dict as a shallow copy.

    '''
    z = x.copy()
    z.update(y)

    if len(z) < len(x) + len(y):
        raise ValueError("Non unique node names!")

    return z
