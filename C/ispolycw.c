def ispolycw(poly):
    """
    Tests if a series of points around the edges of a polygon is wound counterclockwise.
    Returns True for a clockwise polygon, False for a counterclockwise polygon.
    (when looking up at the sky).  This applies to x,y map projection coordinates - if used on
    a polygon in az,el coordinates, the definitions of clockwise and counterclockwise are reversed.
    """
    tot=sum([cross(*vecs) for vecs in zip(poly,poly[1:])])

return tot>0
