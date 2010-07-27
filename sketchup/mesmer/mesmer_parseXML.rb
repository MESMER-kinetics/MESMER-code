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
load Pathname.new(File.dirname(__FILE__)) + 'sphere.rb'
load Pathname.new(File.dirname(__FILE__)) + 'element.rb'
#--------------------

class MesmerObject

  # Instance variables
  # attr_accessor :simulations

  # Initiation procedure
  def initialize
    @simulations = ReactionScheme.new
    @angstrom = 0.1
    @lowestEnergy = 9e23
    @highestEnergy = -9e23
    @tempML = []
    @systemWidth = 0.0
    @speciesWidth = 0.0
    @halfSpeciesW = 0.0
    @systemHeight = 0.0
    @colorScheme = 0
    @ribbonWidth = 10
    @eleList = AtomicElementList.new
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
    # initialize the element list
    @eleList.load_elements

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
      # get the ID and description
      identification = nodeMol.attributes["id"]
      if ml[identification] != nil
        result = UI.messagebox("Redefinition of molecule " + identification + ". Proceeds?", MB_YESNO)
        if result == 7 # No
          UI.messagebox("SketchUp Mesmer Ruby script will now abort.", MB_OK)
          Process.exit
        end
        # proceeds and overwrites the previous molecule setting
      end

      mol = Molecule.new
      mol.parse(nodeMol)

      # ml is a Hash whose [key,value] are [identification, molecule]
      ml[identification] = mol
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
        @simulations.moleculelist[molName].totalZPE = tempEne
        @simulations.moleculelist[molName].tempMLOrder = itemCount.to_f
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

    # For the y displacement of the surfaces, the current idea is initially
    # allowing user to put y coordinates for each species, but the goal is to
    # drag the surfaces along the y-axis and update the surfaces plot according
    # to the average of all the displacements. For example, if we have a reaction
    # scheme:
    # (1) A + B -> C
    # (2) A + B -> D
    # (3) A + B -> E and
    # (4) E -> F
    # (5) F -> G + H
    # Then we drag reaction (1) along the y-axis by 2 units. It will results
    # in C being dragged along the y-axis together with the transition state if
    # any by 2 units, but A + B only by a 2/3 units. (as there are three copies
    # of the same source)
    # In ruby, the coordinates are saved in an shared address to be access
    # by both the reaction ribbon and the molecule.
    # The displacements will then be saved back to the XML file as molecular
    # data for reference.

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
    # Loop the reaction list again to get the dimensions.
    rl.reactions.each { |rxnName, rxn| # Set the surface coordinates from the species data
      tempEne = ml[rxn.rct1].zpe
      tempEne += ml[rxn.rct2].zpe if rxn.rct2
      setEnergyRange(tempEne)

      # It is possible that rct2 is deficientReactant
      rct = rxn.rct1
      rct = rxn.rct2 if rxn.rct2Type == "deficientReactant"
      setCoordinates(rct, tempEne)

      if rxn.transitionState
        tempEne = ml[rxn.transitionState].zpe
        setEnergyRange(tempEne)
        setCoordinates(rxn.transitionState, tempEne)
      end

      tempEne = ml[rxn.pdt1].zpe
      tempEne += ml[rxn.pdt2].zpe if rxn.pdt2
      setEnergyRange(tempEne)

      pdt = rxn.pdt1
      pdt = rxn.pdt2 if !isInTempML(pdt) && rxn.pdt2
      setCoordinates(pdt, tempEne)
    }
    @systemHeight = @highestEnergy - @lowestEnergy

    # The size should be refrained by the ratio of 1:2.
    # If the height of the system is 1, then the width of the system should be about 1.414 times as long.
    @systemWidth = @systemHeight * 2.0
    @speciesWidth = @systemWidth / (@tempML.size - 1).to_f
    @halfSpeciesW = @speciesWidth / 2.0
    @ribbonWidth = 0.03 * @systemHeight
    @angstrom = @ribbonWidth * 2.0
    #=====================================================

    # Draw the surface one by one using Bezier surfaces
    # call draw molecules functions while drawing the surface
    rl.reactions.each{ |rxnName, rxn|
      # decides color for this reaction
      rc_component = rand
      pc_component = rand
      tc_component = rand
      fmu = Sketchup.active_model.materials.add "UpperSurface"
      fmd = Sketchup.active_model.materials.add "LowerSurface"
      fmu.color = [tc_component, rc_component, pc_component]
      fmd.color = [tc_component*0.8, rc_component*0.8, pc_component*0.8]

      surfEntities = Array.new # One reaction is a group
      surfG = Sketchup.active_model.entities.add_group

      pesR1 = ml[rxn.rct1]
      pesR1 = ml[rxn.rct2] if rxn.rct2 && ml[rxn.rct2].tempMLOrder != -1.0
      r_zpe = pesR1.totalZPE
      r_dsp = @speciesWidth * pesR1.tempMLOrder
      moltxt = rxn.rct1.to_s
      moltxt = moltxt + " + " + rxn.rct2.to_s if rxn.rct2
      if !pesR1.isDrawn
        surfG.entities.add_text moltxt, [@speciesWidth * pesR1.tempMLOrder, 0.0, r_zpe - @halfSpeciesW]
        mol = drawMolecule(pesR1, @speciesWidth * pesR1.tempMLOrder, 0.0, r_zpe - 2.0 * @halfSpeciesW)
        #surfEntities.push(mol)
        pesR1.isDrawn = true
      end

      # which product is on pes?
      pesP1 = ml[rxn.pdt1]
      pesP1 = ml[rxn.pdt2] if rxn.pdt2 && ml[rxn.pdt2].tempMLOrder != -1.0
      p_zpe = pesP1.totalZPE
      p_dsp = @speciesWidth * pesP1.tempMLOrder
      moltxt = rxn.rct1.to_s
      moltxt = moltxt + " + " + rxn.pdt2.to_s if rxn.pdt2
      if !pesP1.isDrawn
        surfG.entities.add_text rxn.pdt1.to_s, [@speciesWidth * pesP1.tempMLOrder, 0.0, p_zpe - @halfSpeciesW]
        mol = drawMolecule(pesP1, @speciesWidth * pesP1.tempMLOrder, 0.0, p_zpe - 2.0 * @halfSpeciesW)
        #surfEntities.push(mol)
        pesP1.isDrawn = true
      end

      if rxn.transitionState # This does the ribbon connecting directly from the reactant -> TS -> product
        pesTS = ml[rxn.transitionState]
        ts_zpe = pesTS.zpe
        ts_dsp = pesTS.tempMLOrder

        pt1 = [@speciesWidth * pesR1.tempMLOrder, 0.0, r_zpe]
        pt2 = [@speciesWidth * pesTS.tempMLOrder, 0.0, ts_zpe]
        pMesh1 = Bezier.curveSurface(pt1, pt2, @ribbonWidth)
        surfG.entities.add_faces_from_mesh pMesh1, 12, fmu, fmd

        pt1 = [@speciesWidth * pesTS.tempMLOrder, 0.0, ts_zpe]
        pt2 = [@speciesWidth * pesP1.tempMLOrder, 0.0, p_zpe]
        pMesh2 = Bezier.curveSurface(pt1, pt2, @ribbonWidth)
        surfG.entities.add_faces_from_mesh pMesh2, 12, fmu, fmd

        # draw transition state molecule
        surfG.entities.add_text rxn.transitionState.to_s, [@speciesWidth * pesTS.tempMLOrder, 0.0, ts_zpe + @halfSpeciesW]
        mol = drawMolecule(pesTS, @speciesWidth * pesTS.tempMLOrder, 0.0, ts_zpe + 2.0 * @halfSpeciesW)
        #surfEntities.push(mol)
        pesTS.isDrawn = true

      else # This does the ribbon connecting directly from the reactant to the product
        colorVector = [rc_component, 0.0, pc_component]
        pt1 = [@speciesWidth * pesR1.tempMLOrder, 0.0, r_zpe]
        pt2 = [@speciesWidth * pesP1.tempMLOrder, 0.0, p_zpe]
        pMesh1 = Bezier.curveSurface(pt1, pt2, @ribbonWidth)
        surfG.entities.add_faces_from_mesh pMesh1, 12, fmu, fmd
      end

      #surfEntities.push(surfG)
      #Sketchup.active_model.entities.add_group surfG, surfEntities
    }

    #---------------------------------------
    # Draw a case to be the molecular reservoir
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
    styles = Sketchup.active_model.styles.active_style

  end

  def drawMolecule(mol, ax, ay, az)
    # pass a molecular object in, assuming the coordinates, bonds and atom types
    # are well defined.
    # The plotted molecule should be intact as an object and allow user to rotate it.
    # Todo:

    atoms = mol.atoms
    bonds = mol.bonds

    molGroups = Array.new
    atoms.each { |item|
      molGroups.push(drawAtom(item, ax, ay, az, 2))
    }

    bonds.each {|item|
      molGroups.push(drawBond())
    }

    #return Sketchup.active_model.entities.add_group(molGroups)
  end

  def drawAtom(atom, ax, ay, az,  option)
    # Atom can be a sphere, single dot, or even nothing.
    atomPos = Geom::Point3d.new(atom.x * @angstrom + ax, atom.y * @angstrom + ay, atom.z * @angstrom + az)
    case option
      when 0 # draw nothing
      when 1 # draw a sphere of fixed diameter relative to system size
        # the fixed diameter of the atom is 1 angstrom
        return DrawSphere.create_geometry(atomPos, 1.0 * @angstrom, "Xx")
      when 2 # draw a sphere of van der Waals radius relative to system size
        return DrawSphere.create_geometry(atomPos, Math::sqrt(@eleList[atom.atomName].rcov) * @angstrom, atom.atomName)
    end
  end


  def drawBond()
    # Bond can be single cylinder, double cylinders, sticks, or no lines
    #
    group = nil
    return group
  end

end

