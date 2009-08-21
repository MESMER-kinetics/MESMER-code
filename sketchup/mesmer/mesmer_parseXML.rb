#
#  mesmer_parseXML.rb
#
#
#  Created by Chi-Hsiu Liang on 27/07/2009.
#  Copyright (c) 2009 Chi-Hsiu Liang. All rights reserved.
#

#--------------------
# inclusions specific for this code
require 'rexml/document'
include REXML
require 'pathname' # 1.8
load Pathname.new(File.dirname(__FILE__)) + 'reaction.rb'
load Pathname.new(File.dirname(__FILE__)) + 'bezier.rb'
#--------------------

class MesmerObject

  # Instance variables
  # attr_accessor :simulations

  # Initiation procedure
  def initialize
    @simulations = ReactionScheme.new
    @molecule_scale_factor = 1
    @vanderwaals_scale_factor = 1
    @x_displacement = (100.0).cm
    @lowestEnergy = 9e23
    @highestEnergy = -9e23
    @tempML = []
    @systemWidth = 0.0
    @speciesWidth = 0.0
    @halfSpeciesW = 0.0
    @systemHeight = 0.0
    @colorScheme = 0
    @ribbonWidth = 10
  end

  # Member function which loads the file and calls the parsing procedure
  def load_file(fullFilePath)
    file = File.new(fullFilePath)
    doc = Document.new(file)

    # The next line if uncommented allows the loaded file to be shown on Ruby Console
    #puts doc

    self.parseXMLFile(doc)
    self.drawSurface

    return 0
  end

  # This routine should extract all data from the XML file into a single object
  # for later processes.
  def parseXMLFile(doc)
    # parse the Title section
    @simulations.title = XPath.first( doc, "//title" )
    if ! @simulations.title
      puts "No title."
    else
      puts @simulations.title.to_s
    end

    # create a new instance of moleculeList object
    ml = MoleculeList.new
    doc.elements.each("*/moleculeList/molecule"){ |nodeMol|
      mol = Molecule.new

      # get the ID and description
      identification = nodeMol.attributes["id"]
      if ml.molecules[identification] != nil
        result = UI.messagebox("Redefinition of molecule " +
        identification + ". Proceeds?", MB_YESNO)
        if result == 7 # No
          UI.messagebox("SketchUp Mesmer Ruby script will now abort.", MB_OK)
          Process.exit
        end
        # proceeds and overwrites the previous molecule setting
      end

      mol.description = nodeMol.attributes["description"]

      #-----------------------------------
      # Parse the data inside propertyList
      nodeMol.elements.each("*/property"){ |nodePpt|
        val = nil
        dictRef = nodePpt.attributes["dictRef"]
        text = nodePpt.elements[1].get_text
        val = text.value if text
        if val != nil
          if dictRef == "me:ZPE"
            mol.zpe = val.to_f
          elsif dictRef == "me:rotConsts"
            # replace each tab, carrige return and newline with a space,
            # and then use space as the delimiter to split the string into array
            mol.rotConsts = val.gsub(/[\t\r\n]/,' ').split(" ")
          elsif dictRef == "me:symmetryNumber"
            mol.symm = val.to_f
            #puts "symmetryNumber = " + mol.symm
          elsif dictRef == "me:frequenciesScaleFactor"
            mol.scaleFactor = val.to_f
          elsif dictRef == "me:vibFreqs"
            mol.vibFreq = val.gsub(/[\t\r\n]/,' ').split(" ")
          elsif dictRef == "me:MW"
            mol.molecularWeight = val.to_f
          elsif dictRef == "me:spinMultiplicity"
            mol.spinMultiplicity = val.to_f
          elsif dictRef == "me:electronicExcitation"
            mol.eleExc = val.gsub(/[\t\r\n]/,' ').split(" ")
          end
        end
      }
      #-----------------------------------

      # ml.molecules is a Hash whose [key,value] are [identification, molecule]
      ml.molecules[identification] = mol
    }

    # parse the reactionList
    # create a new instance of moleculeList object
    rl = ReactionList.new
    doc.elements.each("*/reactionList/reaction"){ |nodeRxn|
      rxn = Reaction.new

      # get the ID and description
      identification = nodeRxn.attributes["id"]
      if rl.reactions[identification] != nil
        result = UI.messagebox("Redefinition of reaction " +
        identification + ". Proceeds?", MB_YESNO)
        if result == 7 # No
          UI.messagebox("SketchUp Mesmer Ruby script will now abort.", MB_OK)
          Process.exit
        end
        # proceeds and overwrites the previous reaction setting
      end

      # parse the 1st reactant
      nodeRct1 = nodeRxn.elements[1, 'reactant']
      if nodeRct1
        rxn.rct1 = XPath.first(nodeRct1, "./molecule[1]").attributes["ref"]
        rxn.rct1Type = XPath.first(nodeRct1, "./molecule[1]").attributes["me:type"]
      end

      # parse the 2nd reactant if any
      nodeRct2 = nodeRxn.elements[2, 'reactant']
      if nodeRct2 != nil
        rxn.rct2 = XPath.first(nodeRct2, "./molecule[1]").attributes["ref"]
        rxn.rct2Type = XPath.first(nodeRct2, "./molecule[1]").attributes["me:type"]
      end

      # parse the 1st product if any
      nodePdt1 = nodeRxn.elements[1, 'product']
      if nodePdt1 != nil
        rxn.pdt1 = XPath.first(nodePdt1, "./molecule[1]").attributes["ref"]
        rxn.pdt1Type = XPath.first(nodePdt1, "./molecule[1]").attributes["me:type"]
      end

      # parse the 2nd product if any
      nodePdt2 = nodeRxn.elements[2, 'product']
      if nodePdt2 != nil
        rxn.pdt2 = XPath.first(nodePdt2, "./molecule[1]").attributes["ref"]
        rxn.pdt2Type = XPath.first(nodePdt2, "./molecule[1]").attributes["me:type"]
      end

      # parse the transition state
      nodeTS = nodeRxn.elements[1, 'me:transitionState']
      if nodeTS != nil
        rxn.transitionState = XPath.first(nodeTS, "./molecule[1]").attributes["ref"]
      end

      # check the micro rate calculator
      nodeTxt = nodeRxn.elements[1, 'me:MCRCMethod']
      val = nodeTxt.text if nodeTxt
      rxn.microRateCalculator = val if val
      rl.reactions[identification] = rxn
    }

    # parse the species profile and Bartis-Widom rates.


    # pass the objects back
    @simulations.moleculelist = ml
    @simulations.reactionlist = rl
    return 0
  end

  def setEnergyRange(zpe)
    @lowestEnergy = zpe if zpe < @lowestEnergy
    @highestEnergy = zpe if zpe > @highestEnergy
  end

  # Set the totalZPE and tempMLOrder components of the molecule
  def setCoordinates(molName, tempEne)
    itemCount = 0
    @tempML.each{|item|
      if item == molName
        @simulations.moleculelist.molecules[molName].totalZPE = tempEne
        @simulations.moleculelist.molecules[molName].tempMLOrder = itemCount.to_f
        return 0
      end
      itemCount += 1
    }
    return 1
  end

  def isInTempML(molName)
    @tempML.each{|item|
      if item == molName
        return true
      end
    }
    return false
  end

  # visualize the data collected
  def drawSurface
    # The plan to give some user interactions in sketchup. The interactions
    # would include dragging selected molecules to a existing surface and
    # then modify the reactionlist. So, basically there should be
    # A. a function to write the reaction data back to a XML file (writeXML),
    # B. a function to visualize the reactionlist (drawSurface),
    # C. a function to read the XML file (parseXMLFile),
    # D. a set of functions to operate on the screen objects and manipulate
    #    the data, including drag and drop of molecules, writing interactive
    #    comments on the GUI/XML.

    # Erase all data that has been drawn if this function is called to redraw
    # the surface.
    model_list = Sketchup.active_model.entities
    Sketchup.active_model.entities.each{ |item|
      Sketchup.active_model.entities.erase_entities(item)
    }
    #---------------------------------------
    # Draw surface according to reactionlist
    #---------------------------------------

    rl = @simulations.reactionlist
    ml = @simulations.moleculelist

    # The basic idea of PES sequences is to create an Array of molecular names.
    # In each reaction, molecule names are inserted into the Array according
    # to the reaction specification. For source terms, only the deficient reactant
    # is considered.
    hasTwoRcts = false

    sizeFlag = 0
    parseCount = rl.reactions.size  # the number of reactions remain unparsed
    #=====================================================================
    while parseCount > 0 do

      rl.reactions.each { |rxnName, rxn|

        next if rxn.parsed == true

        #---------------------------------------------
        # parse reactant

        rct = nil
        if rxn.rct2
          #puts "There are two reactants"
          hasTwoRcts = true
          if rxn.rct1Type == "deficientReactant"
            rct = rxn.rct1
          else
            rct = rxn.rct2
          end
        else
          rct = rxn.rct1
        end

        # see if this molecule is in @tempML.
        rctInTempML = false
        rctLoc = 0
        @tempML.each{|item|
          if item == rct
            rctInTempML = true
            break
          end
          rctLoc += 1
        }

        if !rctInTempML
          if hasTwoRcts
            rctLoc = 0
            @tempML.insert(0, rct) # insert the rct to the beginning
            rctInTempML = true
          else
            # if reactant not in the list, need to see if the product
            # is already in the list
          end
        end
        #---------------------------------------------

        #---------------------------------------------
        # parse transition state
        if rxn.transitionState
          @tempML.each{|i1|
            if i1 == rxn.transitionState
              result = UI.messagebox("Error: Reuse of transition state " +
              rxn.transitionState + ". Check input file and correct it.", MB_OK)
              Process.exit
            end
          }
        end
        #---------------------------------------------

        #---------------------------------------------
        # parse product
        pdt1InTempML = pdt2InTempML = false
        pdtLoc = 0
        @tempML.each{ |i2|
          # The difficult situation here is that different reactions can have
          # the same product. For example, OH, can appear in a reaction
          # as a sink together with different counterparts. And, it is possible
          # for two reactions to have exactly the same products but in different
          # order in the XML file.
          if i2 == rxn.pdt1
            pdt1InTempML = true
            break
          elsif i2 == rxn.pdt2
            pdt2InTempML = true
          end
          pdtLoc += 1
        }

        #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx# if both products are not in @tempML
        if !pdt1InTempML && !pdt2InTempML
          if rctInTempML
            # put the TS and product right after the reactant.
            # if it is a sink, push it back into the end of the Array
            if rxn.pdt1Type == "sink"
              @tempML.insert(rctLoc+1, rxn.transitionState) if rxn.transitionState
              @tempML.push(rxn.pdt1)
            else
              if rxn.transitionState
                @tempML.insert(rctLoc+1, rxn.transitionState, rxn.pdt1)
              else
                @tempML.insert(rctLoc+1, rxn.pdt1)
              end
            end
          else
            # if both reactant and product are not in @tempML or if no element
            # is added in final two cycles, put all species in the end of the PES,
            # parse the reaction later
            if parseCount == 1 || sizeFlag > 2
              # in case that this reaction is not related with all other reactions
              # put all species in the end of the PES
              @tempML.push(rxn.rct1)
              @tempML.push(rxn.transitionState) if rxn.transitionState
              @tempML.push(rxn.pdt1)
              sizeFlag = 0
            else
              sizeFlag += 1
              # break the each loop without removing the entry
              break
            end
          end
        #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx# one of the product is in @tempML
        else
          if rctInTempML # if both reactant and one of the products are in @tempML
            tsLoc = (rctLoc + pdtLoc) / 2 + 1
            @tempML.insert(tsLoc, rxn.transitionState) if rxn.transitionState
          else
            @tempML.insert(pdtLoc, rxn.transitionState) if rxn.transitionState
            @tempML.insert(pdtLoc, rct)
          end
        end
        #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        #---------------------------------------------

        # before exiting, remove this entry
        rxn.parsed = true
        parseCount -= 1
        #puts "Marked " + rxnName + " as parsed."
        #puts parseCount.to_s + " entries remain unparsed."
        sizeFlag = 0 # this is to show that the size of rl.reactions is reduced.
      }
    end
    #=====================================================================

    #=====================================================
    puts "The number of species in PES is " + @tempML.size.to_s
    puts "Arranging the coordinates..."

    # Assign horizontal locations to species.
    # The size should be refrained by the ratio of 1:sqrt(2), which is the aspect ratio of an A4 page.
    # If the height of the system is 1, then the width of the system should be about 1.414 times as long.

    # Loop the reaction list again to get the dimensions.
    rl.reactions.each { |rxnName, rxn| # Set the surface coordinates from the species data
      tempEne = ml.molecules[rxn.rct1].zpe
      tempEne += ml.molecules[rxn.rct2].zpe if rxn.rct2
      setEnergyRange(tempEne)

      # It is possible that rct2 is deficientReactant
      rct = rxn.rct1
      rct = rxn.rct2 if rxn.rct2Type == "deficientReactant"
      setCoordinates(rct, tempEne)

      if rxn.transitionState
        tempEne = ml.molecules[rxn.transitionState].zpe
        setEnergyRange(tempEne)
        setCoordinates(rxn.transitionState, tempEne)
      end

      tempEne = ml.molecules[rxn.pdt1].zpe
      tempEne += ml.molecules[rxn.pdt2].zpe if rxn.pdt2
      setEnergyRange(tempEne)

      pdt = rxn.pdt1
      pdt = rxn.pdt2 if !isInTempML(pdt) && rxn.pdt2
      setCoordinates(pdt, tempEne)
    }
    @systemHeight = @highestEnergy - @lowestEnergy
    # Need a line here to do conversion between Hatree and kJ/mol
    @systemWidth = @systemHeight * Math.sqrt(2)
    @speciesWidth = @systemWidth / (@tempML.size - 1).to_f
    @halfSpeciesW = @speciesWidth / 2.0
    @ribbonWidth = 0.03 * @systemHeight
    #=====================================================

    # Draw the surface one by one using Bezier surfaces
    # call draw molecules functions while drawing the surface
    rl.reactions.each{ |rxnName, rxn|
      entities = Sketchup.active_model.entities
      surfaceGroup = Sketchup.active_model.entities.add_group # One reaction is a group

      pesR1 = ml.molecules[rxn.rct1]
      pesR1 = ml.molecules[rxn.rct2] if rxn.rct2 && ml.molecules[rxn.rct2].tempMLOrder != -1.0
      r_zpe = pesR1.totalZPE
      r_dsp = @speciesWidth * pesR1.tempMLOrder
      rc_component = (r_zpe - @lowestEnergy) / @systemHeight
      entities.add_text rxn.rct1.to_s, [@speciesWidth * pesR1.tempMLOrder, 0.0, r_zpe - @halfSpeciesW]

      # which product is on pes?
      pesP1 = ml.molecules[rxn.pdt1]
      pesP1 = ml.molecules[rxn.pdt2] if rxn.pdt2 && ml.molecules[rxn.pdt2].tempMLOrder != -1.0
      p_zpe = pesP1.totalZPE
      p_dsp = @speciesWidth * pesP1.tempMLOrder
      pc_component = (p_zpe - @lowestEnergy) / @systemHeight
      entities.add_text rxn.pdt1.to_s, [@speciesWidth * pesP1.tempMLOrder, 0.0, p_zpe - @halfSpeciesW]

      if rxn.transitionState
        pesTS = ml.molecules[rxn.transitionState]
        ts_zpe = pesTS.zpe
        ts_dsp = pesTS.tempMLOrder
        tc_component = (ts_zpe - @lowestEnergy) / @systemHeight

        f_material = Sketchup.active_model.materials.add "Surface"
        f_material.color = [tc_component, rc_component, pc_component]

        #puts "Color is:" + colorVector.to_s
        pt1 = [@speciesWidth * pesR1.tempMLOrder, 0.0, r_zpe]
        pt2 = [@speciesWidth * pesTS.tempMLOrder, 0.0, ts_zpe]
        pMesh1 = Bezier.curveSurface(pt1, pt2, @ribbonWidth)
        surfaceGroup.entities.add_faces_from_mesh pMesh1, 12, f_material

        entities.add_text rxn.transitionState.to_s, [@speciesWidth * pesTS.tempMLOrder, 0.0, ts_zpe + @halfSpeciesW]

        pt1 = [@speciesWidth * pesTS.tempMLOrder, 0.0, ts_zpe]
        pt2 = [@speciesWidth * pesP1.tempMLOrder, 0.0, p_zpe]
        pMesh2 = Bezier.curveSurface(pt1, pt2, @ribbonWidth)
        surfaceGroup.entities.add_faces_from_mesh pMesh2, 12, f_material
      else
        colorVector = [rc_component, 0.0, pc_component]
        pt1 = [@speciesWidth * pesR1.tempMLOrder, 0.0, r_zpe]
        pt2 = [@speciesWidth * pesP1.tempMLOrder, 0.0, p_zpe]
        pMesh1 = Bezier.curveSurface(pt1, pt2, @ribbonWidth)
        surfaceGroup.entities.add_faces_from_mesh pMesh1, 12, f_material
      end
    }

    #---------------------------------------
    # Draw a case to be the molecular reservoir
    #---------------------------------------

    #---------------------------------------
    # Draw molecules on surface
    #---------------------------------------

    #---------------------------------------
    # Draw molecules on reservior
    #---------------------------------------


    #---------------------------------------
    # Draw comments
    #---------------------------------------

    #---------------------------------------
    # Set scenes for visualization
    #---------------------------------------

    # Set the ortho scene and zoom to extent.
    # This step should be done in the end of all drawing.
    Sketchup.send_action("viewFront:")
    Sketchup.send_action("viewZoomExtents:")
    model = Sketchup.active_model
    view = model.active_view
    camera = view.camera
    status = camera.perspective?
    camera.perspective = false if status
  end

end

