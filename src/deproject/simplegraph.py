"""
A simple/basic graph structure designed to help group together annuli.

:Copyright: Smithsonian Astrophysical Observatory (2018)
:Author: Douglas Burke (dburke@cfa.harvard.edu)
"""

__all__ = ("SimpleGraph", )


class SimpleGraph:
    """A graph where links are bi-directional

    The identifiers must be hashable (as they are used in sets).
    This is not optimised, and so is intended only for small graphs.
    It is expected that a link is manually created for each identifier,
    that is::

        >>> graph = SimpleGraph()
        >>> for n in ['n1', 'n2', 'n3, 'n4', 'n5']:
        ...     graph.add_link(n, n)
        ...
        >>> graph.add_link('n2', 'n3')
        >>> graph.add_link('n3', 'n4')
        >>> graph.get_groups()
        [['n1'], ['n2', 'n3', 'n4'], ['n5']]

    """

    def __init__(self):
        self.groups = []

    def add_link(self, id1, id2):
        """Adds a link between two identifiers (which can be the same).

        The order between id1 and id2 is assumed to be bi-directional.

        Parameters
        ----------
        id1, id2
            identifier of a pair-wise link

        """

        g1 = None
        g2 = None
        for i, grp in enumerate(self.groups):
            if id1 in grp:
                assert g1 is None
                g1 = i

            if id2 in grp:
                assert g2 is None
                g2 = i

        # Not known, so add in. This is not expected to happen,
        # since the idea is to explicitly call with id1 == id2,
        # but support it.
        #
        if g1 is None and g2 is None:
            self.groups.append(set([id1, id2]))
            return

        # Already known about
        if g1 == g2:
            return

        # It is expected that g1 and g2 are both set, but this
        # is not guaranteed.
        #
        if g2 is None:
            self.groups[g1].add(id2)
            return

        if g1 is None:
            self.groups[g2].add(id1)
            return

        # id1 and id2 are part of the g1 and g2 groups, by construction,
        # so they do not need to be added.
        #
        newgrp = self.groups[g1].union(self.groups[g2])
        self.groups[g1] = newgrp
        del self.groups[g2]

    def get_groups(self):
        """Return a copy of the groups, that is linked identifiers.

        Returns
        -------
        groups : list of list
            Each group is identified by the nodes in the list. There is
            no guarantee of ordering of either list.
        """

        out = []
        for grp in self.groups:
            out.append(list(grp))

        return out
