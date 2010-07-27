#
#  mesmer_menuItem.rb
#
#
#  Created by Chi-Hsiu Liang on 27/07/2009.
#  Copyright (c) 2009 Chi-Hsiu Liang. All rights reserved.
#


#--------------------
# reads in the target file specified by the user. Will provide a function (button) to refresh
# the file if it is modified externally. The basic idea is: if a file is loaded, Sketchup
# will open it in a new window and draw the XML file data out in that window. This may
# force other windows to close if run in MS Windows, but shall not be a problem in Mac.
#--------------------
class XMLImporter < Sketchup::Importer

  # This method is called by SketchUp to determine the description that
  # appears in the File > Import dialog's pulldown list of valid
  # importers.
  def description
    return "MESMER XML Importer (*.xml)"
  end

  # This method is called by SketchUp to determine what file extension
  # is associated with your importer.
  def file_extension
    return "xml"
  end

  # This method is called by SketchUp to get a unique importer id.
  def id
    return "com.sketchup.importers.mesmer_xml"
  end

  # This method is called by SketchUp to determine if the "Options"
  # button inside the File > Import dialog should be enabled while your
  # importer is selected.
  def supports_options?
    return true
  end

  # This method is called by SketchUp when the user clicks on the
  # "Options" button inside the File > Import dialog. You can use it to
  # gather and store settings for your importer.
  def do_options
    # In a real use you would probably store this information in an
    # instance variable.
         my_settings = UI.inputbox(['My Import Option:'], ['1'],
           "Import Options")
  end

  # This method is called by SketchUp after the user has selected a file
  # to import. This is where you do the real work of opening and
  # processing the file.
  def load_file(file_path, status)
    n = MesmerObject.new
    n.load_file(file_path)
    return 0 # 0 is the code for a successful import
  end
end

Sketchup.register_importer(XMLImporter.new)