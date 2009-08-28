# Copyright 2004-2005, @Last Software, Inc.

# This software is provided as an example of using the Ruby interface
# to SketchUp.

# Permission to use, copy, modify, and distribute this software for
# any purpose and without fee is hereby granted, provided that the above
# copyright notice appear in all copies.

# THIS SOFTWARE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

module Bezier

  # Evaluate a Bezier curve at a parameter.
  # The curve is defined by an array of its control points.
  # The parameter ranges from 0 to 1
  # This is based on the technique described in "CAGD  A Practical Guide, 4th Editoin"
  # by Gerald Farin. page 60

  def Bezier.eval(pts, t)

    degree = pts.length - 1
    if degree < 1
        return nil
    end

    t1 = 1.0 - t
    fact = 1.0
    n_choose_i = 1

    x = pts[0].x * t1
    y = pts[0].y * t1
    z = pts[0].z * t1

    for i in 1...degree
    fact = fact*t
    n_choose_i = n_choose_i*(degree-i+1)/i
        fn = fact * n_choose_i
    x = (x + fn*pts[i].x) * t1
    y = (y + fn*pts[i].y) * t1
    z = (z + fn*pts[i].z) * t1
    end

    x = x + fact*t*pts[degree].x
    y = y + fact*t*pts[degree].y
    z = z + fact*t*pts[degree].z

    Geom::Point3d.new(x, y, z)

  end # method eval

  # Evaluate the curve at a number of points and return the points in an array
  def Bezier.points(pts, numpts)

    curvepts = []
    dt = 1.0 / numpts

    # evaluate the points on the curve
    for i in 0..numpts
      t = i * dt
      curvepts[i] = Bezier.eval(pts, t)
    end

    curvepts
  end

  # Create a Bezier curve in SketchUp
  def Bezier.curve(pts, numseg = 16)

    model = Sketchup.active_model
    entities = model.active_entities
    model.start_operation "Bezier Curve"

    curvepts = Bezier.points(pts, numseg)
    curvepts
  end

  # This code creates a Bezier surface by creating two parallel Bezier curves
  # and then stitching them together into a surface.
  def Bezier.curveSurface(point1, point2, surfaceWidth) # giving the first point, end point and surface width
    middle = (point1.x + point2.x) / 2.0
    pts1 = []
    pts1.push(point1)
    pts1.push([middle, point1.y, point1.z])
    pts1.push([middle, point2.y, point2.z])
    pts1.push(point2)
    pts2 = []
    pts2.push([point1.x, point1.y + surfaceWidth, point1.z])
    pts2.push([middle, point1.y + surfaceWidth, point1.z])
    pts2.push([middle, point2.y + surfaceWidth, point2.z])
    pts2.push([point2.x, point2.y + surfaceWidth, point2.z])

    #Sketchup.active_model.start_operation "Create surface"
    curve1 = Bezier.curve(pts1, 20)
    curve2 = Bezier.curve(pts2, 20)

    group = Sketchup.active_model.entities.add_group
    pMesh = Geom::PolygonMesh.new
    0.upto(curve1.size - 2) do |a|
      pMesh.add_polygon(curve1[a], curve1[a+1], curve2[a+1],curve2[a])
    end
    return pMesh
  end

end # module Bezier

