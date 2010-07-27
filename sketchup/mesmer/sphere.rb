#
#  sphere.rb
#
#
#  Created by Chi-Hsiu Liang on 27/07/2009.
#  Copyright (c) 2009 Chi-Hsiu Liang. All rights reserved.
#

#
# Ripped and modified based on spheretool.rb 2006 jim.foltz@gmail.com
#

# Permission to use, copy, modify, and distribute this software for
# any purpose and without fee is hereby granted, provided that the above
# copyright notice appear in all copies.

# THIS SOFTWARE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#-----------------------------------------------------------------------------

require 'sketchup.rb'

module DrawSphere
  # Create new geometry when the user has selected two points.
  def DrawSphere.create_geometry(p1, radius, atomName, numsegs = 20)
    puts "The sphere to draw: " + p1.to_s + " " + radius.to_s + " " + atomName
    group = Sketchup.active_model.entities.add_group
    entities = group.entities
    circle1 = entities.add_circle p1, [0,1,0], radius, numsegs
    circle2 = entities.add_circle p1+[0,0,-2*radius], [0,0,1], 2.0 * radius, numsegs
    face1 = entities.add_face circle1
    face1.followme circle2
    entities.erase_entities circle2
    entities.each{|item|
      item.material = atomName if item.typename == "Face"
      item.back_material = atomName if item.typename == "Face"
    }
    return group
  end

end # class Sphere

