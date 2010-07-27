#
#  mesmer_parseXML.rb
#
#
#  Created by Chi-Hsiu Liang on 27/07/2009.
#  Copyright (c) 2009 Chi-Hsiu Liang. All rights reserved.
#
require 'sketchup.rb'
require 'pathname' # 1.8

class AtomicElement
  attr_accessor :ID, :symb, :rcov, :rbo, :rVdW, :maxBnd, :mass, :elNeg, :red, :green, :blue

  def initialize
    @ID = -1
    @symb = ""
    @rcov = 0.0
    @rbo = 0.0
    @rVdW = 0.0
    @maxBnd = 0
    @mass = 0.0
    @elNeg = 0.0
    @red = 0.0
    @green = 0.0
    @blue = 0.0
  end
end

class AtomicElementList

  def initialize
    @elementlist = Hash.new(nil)
  end

  def [](a)   # access operator. The
    return @elementlist[a]
  end

  # This function loads molecular settings from element.txt and sets essential
  # materials (colors) for different atoms.
  def load_elements
    model = Sketchup.active_model
    materials = model.materials
    materials.purge_unused

#    puts "All materials in model:"
#    materials.each{|item|
#      puts item.name
#    }

    @elementlist.clear # make sure it is clean
    File.open(Pathname.new(File.dirname(__FILE__)) + 'element.txt').each { |line|
      next if line[0] == "#"
      element = AtomicElement.new
      arguments      = line.split(" ")
      element.ID     = arguments[0].to_i
      element.symb   = arguments[1].to_s
      element.rcov   = arguments[2].to_f
      element.rbo    = arguments[3].to_f
      element.rVdW   = arguments[4].to_f
      element.maxBnd = arguments[5].to_i
      element.mass   = arguments[6].to_f
      element.elNeg  = arguments[7].to_f
      element.red    = arguments[8].to_f
      element.green  = arguments[9].to_f
      element.blue   = arguments[10].to_f
      @elementlist[element.ID]   = element
      @elementlist[element.symb] = element
      # make sure this array is searchable by both ID (int) and symbol (string)
      # But in element.txt we need only to set the data once.
      m = materials.add element.symb.to_s
      m.color = [element.red, element.green, element.blue]
#      puts "Material " + element.symb + " added to the model."
    }

#    puts "All materials in model:"
#    materials.each{|item|
#      puts item.name
#    }

    if @elementlist.size < 100
      result = UI.messagebox("Error: failed loading element.txt.", MB_OK)
      Process.exit
    end
  end
end