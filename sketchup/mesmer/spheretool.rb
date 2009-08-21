# spheretool.rb 2006 jim.foltz@gmail.com
#
# Ripped directly from linetool.rb distributed with Google SketchUp
#
# Copyright 2006, Google Inc.

# This software is provided as an example of using the Ruby interface
# to SketchUp.

# Permission to use, copy, modify, and distribute this software for
# any purpose and without fee is hereby granted, provided that the above
# copyright notice appear in all copies.

# THIS SOFTWARE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#-----------------------------------------------------------------------------

require 'sketchup.rb'

#-----------------------------------------------------------------------------

# To create a new tool in Ruby, you must define a new class that implements
# the methods for the events that you want to respond to.  You do not have
# to implement methods for every possible event that a Tool can respond to.

# Once you have defined a tool class, you select that tool by creating an
# instance of it and passsing it to Sketchup.active_model.select_tool

# This implementation of a tool tries to be pretty complete to show all
# of the kinds of things that you can do in a tool.  This makes it a little
# complicated.  You should also look at the TrackMouseTool defined in
# utilities.rb for an example of a simpler tool.

# This example shows the implementation of a simple sphere tool.  This tool
# is similar to the Circle tool in SketchUp

class SphereTool

  # This is the standard Ruby initialize method that is called when you create
  # a new object.
  def initialize
    @numsegs = 24
    @ip1 = nil
    @ip2 = nil
    @xdown = 0
    @ydown = 0
  end

  # The activate method is called by SketchUp when the tool is first selected.
  # it is a good place to put most of your initialization
  def activate
    # The Sketchup::InputPoint class is used to get 3D points from screen
    # positions.  It uses the SketchUp inferencing code.
    # In this tool, we will have two points;
    # one for the center, one to calculate the radius
    @ip1 = Sketchup::InputPoint.new
    @ip2 = Sketchup::InputPoint.new
    @ip = Sketchup::InputPoint.new
    @drawn = false

    # This sets the label for the VCB
    #Sketchup::set_status_text $exStrings.GetString("Length"), SB_VCB_LABEL
    #Sketchup::set_status_text $exStrings.GetString("Radius"), SB_VCB_LABEL
    Sketchup::set_status_text $exStrings.GetString("Segments"), SB_VCB_LABEL
    Sketchup::set_status_text $exStrings.GetString("#{@numsegs}"), SB_VCB_VALUE

    self.reset(nil)
  end

  # deactivate is called when the tool is deactivated because
  # a different tool was selected
  def deactivate(view)
     view.invalidate if @drawn
  end

  # The onMouseMove method is called whenever the user moves the mouse.
  # because it is called so often, it is important to try to make it efficient.
  # In a lot of tools, your main interaction will occur in this method.
  def onMouseMove(flags, x, y, view)
    if( @state == 0 )
      # We are getting the center of the sphere. Call the pick method
      # on the InputPoint to get a 3D position from the 2D screen position
      # that is passed as an argument to this method.
      @ip.pick view, x, y
      if( @ip != @ip1 )
        # if the point has changed from the last one we got, then
        # see if we need to display the point.  We need to display it
        # if it has a display representation or if the previous point
        # was displayed.  The invalidate method on the view is used
        # to tell the view that something has changed so that you need
        # to refresh the view.
        view.invalidate if( @ip.display? or @ip1.display? )
        @ip1.copy! @ip

        # set the tooltip that should be displayed to this point
        view.tooltip = @ip1.tooltip
      end
    else
      # Getting the radius
      # If you pass in another InputPoint on the pick method of InputPoint
      # it uses that second point to do additional inferencing such as
      # parallel to an axis.
      @ip2.pick view, x, y, @ip1
      view.tooltip = @ip2.tooltip if( @ip2.valid? )
      view.invalidate

      # Update the length displayed in the VCB
      if( @ip2.valid? )
        length = @ip1.position.distance(@ip2.position)
        Sketchup::set_status_text length.to_s, SB_VCB_VALUE
      end

      # Check to see if the mouse was moved far enough to create a line.
      # This is used so that you can create a line by either draggin
      # or doing click-move-click
      if( (x-@xdown).abs > 10 || (y-@ydown).abs > 10 )
        @dragging = true
      end
    end
  end

  # The onLButtonDOwn method is called when the user presses the left mouse button.
  def onLButtonDown(flags, x, y, view)
     # When the user clicks the first time, we switch to getting the
     # second point.  When they click a second time we create the sphere.
     if( @state == 0 )
       @ip1.pick view, x, y
       if( @ip1.valid? )
         @state = 1
         #Sketchup::set_status_text $exStrings.GetString("Select second end"), SB_PROMPT
         Sketchup::set_status_text $exStrings.GetString("Type radius, or select radius end point"), SB_PROMPT
         Sketchup::set_status_text $exStrings.GetString("Radius"), SB_VCB_LABEL
         @xdown = x
         @ydown = y
       end
     else
        # create the sphere on the second click
        if( @ip2.valid? )
           self.create_geometry(@ip1.position, @ip2.position,view)
           self.reset(view)
        end
     end

     # Clear any inference lock
     view.lock_inference
  end

  # The onLButtonUp method is called when the user releases the left mouse button.
  def onLButtonUp(flags, x, y, view)
    # If we are doing a drag, then create the line on the mouse up event
    if( @dragging && @ip2.valid? )
      self.create_geometry(@ip1.position, @ip2.position,view)
      self.reset(view)
    end
  end

  # onKeyDown is called when the user presses a key on the keyboard.
  # We are checking it here to see if the user pressed the shift key
  # so that we can do inference locking
  def onKeyDown(key, repeat, flags, view)
    if( key == CONSTRAIN_MODIFIER_KEY && repeat == 1 )
      @shift_down_time = Time.now

      # if we already have an inference lock, then unlock it
      if( view.inference_locked? )
        # calling lock_inference with no arguments actually unlocks
        view.lock_inference
      elsif( @state == 0 && @ip1.valid? )
        view.lock_inference @ip1
      elsif( @state == 1 && @ip2.valid? )
        view.lock_inference @ip2, @ip1
      end
    end
  end

  # onKeyUp is called when the user releases the key
  # We use this to unlock the interence
  # If the user holds down the shift key for more than 1/2 second, then we
  # unlock the inference on the release.  Otherwise, the user presses shift
  # once to lock and a second time to unlock.
  def onKeyUp(key, repeat, flags, view)
    if( key == CONSTRAIN_MODIFIER_KEY &&
      view.inference_locked? &&
      (Time.now - @shift_down_time) > 0.5 )
      view.lock_inference
    end
  end

  # onUserText is called when the user enters something into the VCB
  # In this implementation, we create a line of the entered length if
  # the user types a length while selecting the second point
  def onUserText(text, view)
    # We only accept input when the state is 1 (i.e. getting the second point)
    # This could be enhanced to also modify the last line created if a length
    # is entered after creating a line.
    #return if not @state == 1
    if @state == 0
      @numsegs = text.to_l
      n, r = @numsegs.divmod(2)
      @numsegs += 1 if r != 0
      Sketchup::set_status_text "", SB_VCB_VALUE
      Sketchup::set_status_text $exStrings.GetString("Radius"), SB_VCB_LABEL
    end
    return if not @ip2.valid?

    # The user may type in something that we can't parse as a length
    # so we set up some exception handling to trap that
    begin
      value = text.to_l
    rescue
      # Error parsing the text
      UI.beep
      puts "Cannot convert #{text} to a Length"
      value = nil
      Sketchup::set_status_text "", SB_VCB_VALUE
    end
    return if !value

    # Compute the direction and the second point
    pt1 = @ip1.position
    vec = @ip2.position - pt1
    if( vec.length == 0.0 )
      UI.beep
      return
    end
    vec.length = value
    pt2 = pt1 + vec

    # Create a line
    self.create_geometry(pt1, pt2, view)
    self.reset(view)
  end

  # The draw method is called whenever the view is refreshed.  It lets the
  # tool draw any temporary geometry that it needs to.
  def draw(view)
    if( @ip1.valid? )
      if( @ip1.display? )
        @ip1.draw(view)
        @drawn = true
      end

      if( @ip2.valid? )
        @ip2.draw(view) if( @ip2.display? )

        # The set_color_from_line method determines what color
        # to use to draw a line based on its direction.  For example
        # red, green or blue.
        view.set_color_from_line(@ip1, @ip2)
        self.draw_geometry(@ip1.position, @ip2.position, view)
        @drawn = true
      end
    end
  end

  # onCancel is called when the user hits the escape key
  def onCancel(flag, view)
    self.reset(view)
  end


  # The following methods are not directly called from SketchUp.  They are
  # internal methods that are used to support the other methods in this class.

  # Reset the tool back to its initial state
  def reset(view)
    # This variable keeps track of which point we are currently getting
    @state = 0

    # Display a prompt on the status bar
    Sketchup::set_status_text($exStrings.GetString("Enter segments [#{@numsegs}], or select center point."), SB_PROMPT)
    #Sketchup::set_status_text $exStrings.GetString("Segments"), SB_VCB_LABEL

    # clear the InputPoints
    @ip1.clear
    @ip2.clear

    if( view )
      view.tooltip = nil
      view.invalidate if @drawn
    end

    @drawn = false
    @dragging = false
  end

  # Create new geometry when the user has selected two points.
  def create_geometry(p1, p2, view)
     model = view.model
     model.start_operation("Draw Sphere")
     radius = p2.distance(p1)
     group = view.model.entities.add_group
     entities = group.entities
     #entities = model.entities
     circle1 = entities.add_circle p1, [0,1,0], radius, @numsegs
     circle2 = entities.add_circle p1+[0,0,-2*radius], [0,0,1], 2.0 * radius, @numsegs
     face1 = entities.add_face circle1
     #face2 = entities.add_face circle2
     face1.followme circle2
     entities.erase_entities circle2
     #view.model.entities.add_cline(p1,p2)
     model.commit_operation
  end

=begin
   def create_geometry(p1, p2, view)
      puts "cg"
      model = view.model
      group = model.entities.add_group
      entities = group.entities
      radius = p2.distance(p1)
      dtheta = dphi = 360.0  / @numsegs
      -90.step(90, dtheta) do |theta|
         0.step(360, dphi) do |phi|
            pts = []
            x = radius * Math::cos(theta.degrees) * Math::cos(phi.degrees)
            y =radius* Math::cos(theta.degrees) * Math::sin(phi.degrees)
            z = radius*Math::sin(theta.degrees)
            pts << (p1 + [x, y, z])
            x = radius*Math::cos((theta+dtheta).degrees) * Math::cos(phi.degrees)
            y = radius*Math::cos((theta+dtheta).degrees) * Math::sin(phi.degrees)
            z = radius*Math::sin((theta+dtheta).degrees)
            pts << (p1 + [x, y, z])
            x = radius*Math::cos((theta+dtheta).degrees) * Math::cos((phi+dphi).degrees)
            y = radius*Math::cos((theta+dtheta).degrees) * Math::sin((phi+dphi).degrees)
            z = radius*Math::sin((theta+dtheta).degrees)
            pts << (p1 + [x, y, z])
            if theta > -90 && theta < 90
               x = radius * Math::cos(theta.degrees) * Math::cos((phi+dphi).degrees)
               y = radius * Math::cos(theta.degrees) * Math::sin((phi+dphi).degrees)
               z = radius * Math::sin(theta.degrees)
            pts << (p1 + [x, y, z])
            end
            #pts.uniq!
            begin
            entities.add_face pts
            rescue
            end
         end
      end
   end
=end

  # Draw the geometry
  def draw_geometry(pt1, pt2, view)
    view.draw_line(pt1, pt2)
    #view.draw(GL_LINES,pt1, pt2)
  end

end # class SphereTool
unless file_loaded?("spheretool.rb")
  UI.menu("Draw").add_item("Sphere") { Sketchup.active_model.select_tool SphereTool.new }
  file_loaded("spheretool.rb")
end
