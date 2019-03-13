"""
Limited testing of deproject.simplegraph

:Copyright: Smithsonian Astrophysical Observatory (2019)
:Author: Douglas Burke (dburke@cfa.harvard.edu)
"""

import pytest

from deproject.simplegraph import SimpleGraph


def _sort_groups(gs):
    return sorted([sorted(g) for g in gs])


def test_empty():
    """No nodes."""

    g = SimpleGraph()
    groups = g.get_groups()
    assert groups == []


def test_0d():
    """Everything is one."""

    g = SimpleGraph()
    g.add_link(23, 23)
    groups = g.get_groups()
    assert groups == [[23]]


def test_0d_really():
    """Everything is one even if we like to repeat."""

    g = SimpleGraph()
    g.add_link(' ', ' ')
    g.add_link(' ', ' ')
    groups = g.get_groups()
    assert groups == [[' ']]


def test_single_link_same_type():
    """A single link (same type)"""

    g = SimpleGraph()
    g.add_link('x ', ' y')
    groups = _sort_groups(g.get_groups())
    assert groups == [[' y', 'x ']]


def test_single_link_diff_type():
    """A single link with different types for the nodes"""

    g = SimpleGraph()
    g.add_link(' x ', True)

    # can not sort this list, so explicitly check
    groups = g.get_groups()
    assert len(groups) == 1
    assert len(groups[0]) == 2
    assert ' x ' in groups[0]
    assert True in groups[0]


def test_odd_nodes_unconnected():
    """Self link plus another link, unconnected (and test node types)"""

    g = SimpleGraph()
    g.add_link('', '')
    g.add_link(300, 2)

    # can not sort this list, so explicitly check
    groups = g.get_groups()
    assert len(groups) == 2
    g1 = groups[0]
    g2 = groups[1]

    def pred1(g):
        """Is the group the ''->'' link?"""
        return len(g) == 1 and g[0] == ''

    def pred2(g):
        """Is the group the 2->300 link?"""
        return len(g) == 2 and 2 in g and 300 in g

    assert (pred1(g1) and pred2(g2)) or (pred1(g2) and pred2(g1))


def test_odd_nodes_connected():
    """Connect up nodes of different types"""

    g = SimpleGraph()
    g.add_link(True, 5)
    g.add_link(5, False)

    groups = g.get_groups()
    assert len(groups) == 1
    assert len(groups[0]) == 3
    assert True in groups[0]
    assert False in groups[0]
    assert 5 in groups[0]


# Could parametrize many of the tests since using the same
# basic data, but leave separate for now.
#
def test_no_pairs():
    """no pairs"""

    g = SimpleGraph()

    for i in range(5, 11):
        g.add_link(i, i)

    groups = g.get_groups()
    assert len(groups) == 6

    vs = _sort_groups(groups)
    assert vs == [[5], [6], [7], [8], [9], [10]], vs


@pytest.mark.parametrize("n1,n2", [(8, 9), (9, 8)])
def test_one_pair(n1, n2):
    """one pair"""

    g = SimpleGraph()

    for i in range(5, 11):
        g.add_link(i, i)

    g.add_link(n1, n2)

    groups = g.get_groups()
    assert len(groups) == 5

    vs = _sort_groups(groups)
    assert vs == [[5], [6], [7], [8, 9], [10]], vs


@pytest.mark.parametrize("n1,n2,m1,m2", [(8, 9, 7, 8),
                                         (9, 8, 7, 8),
                                         (8, 9, 8, 7),
                                         (9, 8, 8, 7)])
def test_two_pairs(n1, n2, m1, m2):
    """two pairs"""

    g = SimpleGraph()

    for i in range(5, 11):
        g.add_link(i, i)

    g.add_link(n1, n2)
    g.add_link(m1, m2)

    groups = g.get_groups()
    assert len(groups) == 4

    vs = _sort_groups(groups)
    assert vs == [[5], [6], [7, 8, 9], [10]], vs
