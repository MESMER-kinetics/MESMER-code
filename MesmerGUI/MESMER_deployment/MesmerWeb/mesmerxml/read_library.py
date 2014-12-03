# -*- coding:utf-8 -*-
__author__ = 'ciaosu'

from lxml import etree
from shortcut import *
import json

def read_library(filename):
    with open(filename, "r") as f:
        lib_xml = f.read()
    try:
        root = etree.fromstring(lib_xml)
    except:
        pass

    moleculelist_node = root
    moleculelist = {}
    for molecule_node in moleculelist_node:
        molecule = {}
        if isinstance(molecule_node, etree._Comment):
            continue
        propertylist_node = getNode(molecule_node, "propertyList")
        if propertylist_node is None:
            continue
        for property_node in propertylist_node:
            if isinstance(molecule_node, etree._Comment):
                continue
            try:
                ## get key ##
                key = getAttr(property_node, "dictRef")
                if key:
                    key = key.replace("me:", "")
                else:
                    key = getAttr(property_node, "title")

                ## get value ##
                data = {}
                value_node = property_node[0]
                tag = etree.QName(value_node).localname
                units = getAttr(value_node, "units")
                value = []
                for i in value_node.text.split(" "):
                    try:
                        value.append(float(i))
                    except:
                        pass
                if units is not None:
                    data["units"] = units
                else:
                    data["units"] = ""
                if tag == "array":
                    data["value"] = value
                else:
                    data["value"] = value[0]
                molecule[key] = data
            except:
                pass
        DOSCMethod_node = getNode(molecule_node, "DOSCMethod")
        if DOSCMethod_node is not None:
            molecule["DOSCMethod"] = DOSCMethod_node.text
        molecule["id"] = getAttr(molecule_node, "id")
        moleculelist[molecule["id"]] = molecule

    with open("molecule_library.json", "w") as f:
        f.write(json.dumps(moleculelist))



read_library("H:\mesmer\librarymols.xml")
