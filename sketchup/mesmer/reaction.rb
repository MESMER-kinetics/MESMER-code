#
#  reaction.rb
#
#
#  Created by Chi-Hsiu Liang on 27/07/2009.
#  Copyright (c) 2009 Chi-Hsiu Liang. All rights reserved.
#
require 'rexml/document'
include REXML
require 'pathname' # 1.8
load Pathname.new(File.dirname(__FILE__)) + 'molecule.rb'

#----------------
class Reaction
  #These are molecule IDs
  attr_accessor :rct1, :rct2, :pdt1, :pdt2, :parsed,
                :rct1Type, :rct2Type, :pdt1Type, :pdt2Type,
                :transitionState, :microRateCalculator


  def initialize
    @rct1 = nil
    @rct2 = nil
    @pdt1 = nil
    @pdt2 = nil
    @parsed = false
    @rct1Type = nil
    @rct2Type = nil
    @pdt1Type = nil
    @pdt2Type = nil
    @transitionState = nil
    @microRateCalculator = nil
  end

  def parse(nodeRxn)
    # pass the node in and populate the member data.
  end
end

#----------------
class ReactionList
  attr_accessor :reactions,:activeCount

  def initialize
    @reactions = Hash.new(nil)
    @activeCount = nil
  end
end

#----------------
class PTs
  attr_accessor :pres, :temp, :unit

  def initialize(pres = 0.0, temp = 0.0, unit = "Torr")
    @pres = pres
    @temp = temp
    @unit = unit
  end

end

#----------------
class Condition
  attr_accessor :bathGas, :pts

  def initialize(bathGas = nil)
    @bathGas = bathGas
    @pts = PTs.new
  end

  def add_pt(singlePoint)
    pts.push(singlePoint)
  end
end

#----------------
class ReactionScheme
  attr_accessor :title, :cond, :reactionlist, :moleculelist #, :modelparameter, :control

  def initialize
    @title = nil
    @cond = Condition.new
    @reactionlist = ReactionList.new
    @moleculelist = MoleculeList.new
    #@modelparameter = nil
    #@control = nil
  end
end