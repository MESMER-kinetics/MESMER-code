#
#  mesmer_data1.rb
#
#
#  Created by Chi-Hsiu Liang on 27/07/2009.
#  Copyright (c) 2009 Chi-Hsiu Liang. All rights reserved.
#

#-----------------------------------------------------------------------------
# The function of this script is to use SketchUp as an interface to manipulate
# and visualize MESMER XML documents, which includes the following
# functionalities:
# 1. Read and visualize MESMER single or multiple SketchUp files.
#   a. PES
#   b. DOS
#   c. k(E)s in grains and possibly cells
#   d. Bartis-Widom rates and species profile
#   e. grid search and fitting results
#   All above both shown in 2D and 3D, which can be achieved by using
#   SketchUp "Scenes".
# 2. Provide a set of UI's for entering data and save to inidividual XML files.
#    This includes an array of file names and an array of target models.
#    The data can include:
#   a. molecular definitions (atoms, coordinates, vibrational frequencies, etc.)
#   b. reaction schemes
#   c. conditions
#   d. 
#-----------------------------------------------------------------------------

#-------------------------------------------------
# Temporary code to call this script from Sketchup
#   Windows:
# load 'M:\codes\Mesmer\sketchup\mesmer_data.rb'
#   or if the code is place in sketchup plugin folder
# load 'C:\Program Files\Google\Google SketchUp 7\Plugins\mesmer_data.rb'
#
#   Mac:
# load '/Users/chi-hsiuliang/codes/sketchup/mesmer_data/mesmer_data.rb'
#-------------------------------------------------
# Do not change the script order in this file

require 'sketchup.rb'
require 'extensions.rb'
require 'LangHandler.rb'

# switch between Windows and Mac Rubies.
if Object::RUBY_PLATFORM =~ /mswin/i
  # The line behind should be the address the user install their full version of Ruby
  # SketchUp provides its own set of Ruby but it is only a minimal setting.
  # Full set of Ruby comes with rexml by default.
  $: << 'c:\\ruby\\lib\\ruby\\1.8'
elsif Object::RUBY_PLATFORM =~ /darwin/i
  $: << "/System/Library/Frameworks/Ruby.framework/Versions/1.8/usr/lib/ruby/1.8/"
end

require 'pathname' # 1.8
load Pathname.new(File.dirname(__FILE__)) + 'mesmer/mesmer_parseXML.rb'
require Pathname.new(File.dirname(__FILE__)) + 'mesmer/mesmer_menuItem.rb'
