#
#  atom.rb
#
#
#  Created by Chi-Hsiu Liang on 27/07/2009.
#  Copyright (c) 2009 Chi-Hsiu Liang. All rights reserved.
#

#Simple definitions of Molecule related objects.

#----------------
class Atom
  attr_accessor :atomicOrder, :atomName, :x , :y, :z

  def initialize
    @atomicOrder = 0
    @atomName = nil
    @x = 0.0
    @y = 0.0
    @z = 0.0
  end

  def r_ij(atom_two)
    r_squared = (@x-atom_two.x)**2 + (@y-atom_two.y)**2 + (@z-atom_two.z)**2
    return Math.sqrt(r_squared)
  end
end


#----------------
class Bond
  attr_accessor :atom1, :atom2, :bondOrder

  def initialize
    @atom1 = @atom2 = @bondOrder = -1
  end

end

#----------------
class Molecule
  attr_accessor :tempMLOrder, :totalZPE, :depth, :description, :sigma, :epsilon, :zpe,
                :rotConsts, :symm, :zpe, :scaleFactor, :spinMultiplicity,
                :eleExc, :vibFreq, :grainEne, :grainDOS, :imFreq,
                :initPopulation, :eqFraction, :deltaEdownExponent,
                :deltaEdownRefTemp, :deltaEdown, :resLoc, :atoms, :bonds,
                :molecularWeight

  # The variables below, if possible, should be properly initialized with a predefined library
  def initialize
    @tempMLOrder = -1.0 # the order of the molecule in tempML -1 if not on tempML
    @totalZPE = nil
    @depth = 0.0

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
    @atoms = []
    @bonds = []
    @molecularWeight = 0.0
  end
end

#----------------
class MoleculeList
  attr_accessor :molecules

  def initialize
    @molecules = Hash.new(nil)
  end
end