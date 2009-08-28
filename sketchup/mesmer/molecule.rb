#
#  atom.rb
#
#
#  Created by Chi-Hsiu Liang on 27/07/2009.
#  Copyright (c) 2009 Chi-Hsiu Liang. All rights reserved.
#

#Simple definitions of Molecule related objects.
require 'rexml/document'
include REXML

#----------------
class Atom
  attr_accessor :index, :atomName, :x , :y, :z

  def initialize
    @index = ""
    @atomName = nil
    @x = 0.0
    @y = 0.0
    @z = 0.0
  end

  def r_ij(atom_two)
    r_squared = (@x-atom_two.x)**2 + (@y-atom_two.y)**2 + (@z-atom_two.z)**2
    return Math.sqrt(r_squared)
  end

  def parse(nodeAtom)
    # pass the node in and populate the member data.
    @index = nodeAtom.attributes["id"].to_s
    @atomName = nodeAtom.attributes["elementType"].to_s
    @x = nodeAtom.attributes["x3"].to_f
    @y = nodeAtom.attributes["y3"].to_f
    @z = nodeAtom.attributes["z3"].to_f
    puts @x.to_s + " " + @y.to_s + " " + @z.to_s
  end
end


#----------------
class Bond
  attr_accessor :atom1, :atom2, :bondOrder

  def initialize
    @atom1 = ""
    @atom2 = ""
    @bondOrder = -1
  end

  def parse(nodeBond)
    # pass the node in and populate the member data.
    atrefs = nodeBond.attributes["atomRefs2"].to_s.split(" ")
    @atom1 = atrefs[0]
    @atom2 = atrefs[1]
    @bondOrder = nodeBond.attributes["order"].to_i
  end
end

#----------------
class Molecule
  attr_accessor :tempMLOrder, :totalZPE, :depth, :isDrawn, :description, :sigma, :epsilon, :zpe,
                :rotConsts, :symm, :zpe, :scaleFactor, :spinMultiplicity,
                :eleExc, :vibFreq, :grainEne, :grainDOS, :imFreq,
                :initPopulation, :eqFraction, :deltaEdownExponent,
                :deltaEdownRefTemp, :deltaEdown, :resLoc, :structureDefined, :atoms, :bonds,
                :molecularWeight

  # The variables below, if possible, should be properly initialized with a predefined library
  def initialize
    @tempMLOrder = -1.0 # the order of the molecule in tempML -1 if not on tempML
    @totalZPE = nil
    @depth = 0.0
    @isDrawn = false

    # Molecule
    @description = nil

    # Bath
    @sigma = 0.0
    @epsilon = 0.0

    # density of states
    @rotConsts = []
    @symm = 1
    @zpe = 0.0
    @scaleFactor = 1.0
    @spinMultiplicity = 1
    @eleExc = []
    @vibFreq = []
    @grainEne = []
    @grainDOS = []

    # transition states
    @imFreq = 0.0

    # population
    @initPopulation = 0.0
    @eqFraction = 0.0

    # well properties
    @deltaEdownExponent = 0.0
    @deltaEdownRefTemp = 0.0
    @deltaEdown = 0.0
    @resLoc = 0.0

    # structure
    @structureDefined = false
    @atoms = []
    @bonds = []
    @molecularWeight = 0.0
  end

  def parse(nodeMol)
    @description = nodeMol.attributes["description"]

    #-----------------------------------
    # Parse the data inside propertyList
    nodeMol.elements.each("*/property"){ |nodePpt|
      val = nil
      dictRef = nodePpt.attributes["dictRef"]
      text = nodePpt.elements[1].get_text
      val = text.value if text
      if val != nil
        if dictRef == "me:ZPE"
          @zpe = val.to_f
        elsif dictRef == "me:rotConsts"
          # replace each tab, carrige return and newline with a space,
          # and then use space as the delimiter to split the string into array
          @rotConsts = val.gsub(/[\t\r\n]/,' ').split(" ")
        elsif dictRef == "me:symmetryNumber"
          @symm = val.to_f
          #puts "symmetryNumber = " + @symm
        elsif dictRef == "me:frequenciesScaleFactor"
          @scaleFactor = val.to_f
        elsif dictRef == "me:vibFreqs"
          @vibFreq = val.gsub(/[\t\r\n]/,' ').split(" ")
        elsif dictRef == "me:MW"
          @molecularWeight = val.to_f
        elsif dictRef == "me:spinMultiplicity"
          @spinMultiplicity = val.to_f
        elsif dictRef == "me:electronicExcitation"
          @eleExc = val.gsub(/[\t\r\n]/,' ').split(" ")
        end
      end
    }
    #-----------------------------------
    #-----------------------------------
    # Parse the data inside atomArray
    nodeMol.elements.each("*/atom"){ |nodeAtom|
puts "I am here 3"
      atom = Atom.new
      atom.parse(nodeAtom)
      @atoms.push(atom)
    }
    #-----------------------------------

    #-----------------------------------
    # Parse the data inside bondArray
    nodeMol.elements.each("*/bond"){ |nodeBond|
      bond = Bond.new
      bond.parse(nodeBond)
      @bonds.push(bond)
    }
    #-----------------------------------

  end
end

#----------------
class MoleculeList
  attr_accessor :molecules

  def initialize
    @molecules = Hash.new(nil)
  end

  def [](a)   # access operator. The
    return @molecules[a]
  end
  
  def []=(key, value)
    @molecules[key] = value
  end
end