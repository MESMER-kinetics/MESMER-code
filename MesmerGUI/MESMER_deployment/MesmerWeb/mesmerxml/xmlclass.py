# -*- coding: utf-8 -*-
__author__ = 'ciaosu'

from lxml import etree
from shortcut import *

class MesmerXML:
    control_options = [
        'printCellDOS',
        'printCellTransitionStateFlux',
        'printReactionOperatorColumnSums',
        'printGrainBoltzmann',
        'printGrainDOS',
        'printGrainkbE',
        'printGrainkfE',
        'printTSsos',
        'printGrainedSpeciesProfile',
        'printGrainTransitionStateFlux',
        'printReactionOperatorSizematrix',
        'printSpeciesProfile',
        'printTunnelingCoefficients',
        'printTunnellingCoefficients',
        'printCrossingCoefficients',
        'MaximumEvolutionTime',
        'testDOS',
        'testMicroRates',
        'testRateConstants',
        'useTheSameCellNumberForAllConditions',
        'hideInactive'
    ]
    def __init__(self, xml):
        try:
            self.root = etree.fromstring(xml)
        except:
            self.error = "syntax error"

        try:
            self.me_ns = self.root.nsmap['me']
            self.ns = self.root.nsmap[None]
            self.xsi_ns = self.root.nsmap['xsi']
        except:
            pass

        self.molecules = {}

    # read condition from xml to json data
    def read_condition(self, condition_node):
        conditions = {}

        # read bathGas
        try:
            bath_gas_node = getNode(condition_node, 'bathGas')
            conditions['bathGas'] = bath_gas_node.text
            self.molecules[conditions['bathGas']]['type'] = 'bathGas'
        except:
            pass


        # read InitialPopulation
        try:
            init_population_node = getNode(condition_node, 'InitialPopulation')
            conditions['initialPopulation'] = getAttr(init_population_node[0], 'population')
            conditions['initialPopulationRef'] = getAttr(init_population_node[0], 'ref')
        except:
            pass

        # read PTs
        try:
            conditions['PTs'] = []
            pts_node = getNode(condition_node, 'PTs')
            for pt in pts_node:
                if isinstance(pt, etree._Comment):
                    continue
                conditions['PTs'].append({
                    'P': getAttr(pt, 'P'),
                    'T': getAttr(pt, 'T')
                })
                conditions['units'] = getAttr(pt, 'units')
        except:
            pass

        return conditions

    # read control from xml to json data
    def read_control(self, control_node):
        control = {
            'options': {}
        }
        for node in control_node:
            if isinstance(node, etree._Comment):
                continue
            else:
                tag = etree.QName(node).localname
                if tag == 'eigenvalues':
                    control['eigenvalues'] = node.text
                elif tag == 'calcMethod':
                    control['calcMethod'] = getAttr(node, 'name', True)
                else:
                    try:
                        self.control_options.index(tag)
                        control['options'][tag] = True
                    except:
                        pass
        return control

    # read modelParameters from xml to json data
    def read_model_parameters(self, parameter_node):
        model_parameters = {}
        for node in parameter_node:
            if isinstance(node, etree._Comment):
                continue
            else:
                parameter_el_tag = etree.QName(node).localname
                model_parameters[parameter_el_tag] = node.text
        return model_parameters

    # read molecule list from xml to data
    def read_molecule_list(self, molecule_list_node):
        molecule_list = []
        for node in molecule_list_node:
            if isinstance(node, etree._Comment):
                continue
            node_tag = etree.QName(node).localname
            if node_tag == "molecule":
                molecule = {
                    'id': getAttr(node, 'id')
                }
                try:
                    # read properties
                    property_list = getNode(node, 'propertyList')
                    for p in property_list:
                        if isinstance(p, etree._Comment):
                            continue
                        data_node = getNode(p, 'scalar')
                        if data_node is None:
                            data_node = getNode(p, 'array')
                        key = p.get('dictRef').replace('me:', '')
                        value = data_node.text
                        molecule[key] = value
                        if key == 'ZPE':
                            molecule['ZPE_unit'] = getAttr(data_node, 'units')

                    # read DOSCMethod
                    dosc_method_node = getNode(node, 'DOSCMethod')
                    molecule['DOSCMethod'] = getAttr(dosc_method_node, 'name', from_text=True)

                    # read deltaEDown
                    if molecule.get('deltaEDown') is None:
                        deltaEDown_node = getNode(node, 'deltaEDown')
                        if deltaEDown_node is not None:
                            molecule['deltaEDown'] = deltaEDown_node.text

                    # read deltaEDownTExponent
                    if molecule.get('deltaEDownTExponent') is None:
                        deltaEDownTExponent_node = getNode(node, 'deltaEDownTExponent')
                        if deltaEDownTExponent_node is not None:
                            molecule['deltaEDownTExponent'] = deltaEDownTExponent_node.text
                            molecule['referenceTemperature'] = getAttr(deltaEDownTExponent_node, 'referenceTemperature')
                except:
                    pass
                if molecule.get('vibFreqs'):
                    molecule['vibFreqs'] = separate(molecule['vibFreqs'], ' ')
                if molecule.get('rotConsts'):
                    molecule['rotConsts'] = separate(molecule['rotConsts'], ' ')
                molecule_list.append(molecule)
                self.molecules[molecule['id']] = molecule
        return molecule_list


    # read molecule list from xml to data
    def read_reaction_list(self, reaction_node):
        reaction_list = []

        for node in reaction_node:
            if isinstance(node, etree._Comment):
                continue
            reaction = {}
            # read id
            reaction['id'] = getAttr(node, 'id')

            # read reactant
            try:
                reactant_list = getNodes(node, 'reactant')
                for index, value in enumerate(reactant_list):
                    m_type = getAttr(value[0], 'type')
                    if m_type is None:
                        m_type = getAttr(value[0], 'role')
                    m_ref = getAttr(value[0], 'ref')
                    reaction['R' + str(index+1) + 'Type'] = m_type
                    reaction['R' + str(index+1) + 'Ref'] = m_ref
                    self.molecules[m_ref]['type'] = m_type
            except:
                pass

            # read product
            try:
                product_list = getNodes(node, 'product')
                for index, value in enumerate(product_list):
                    m_type = getAttr(value[0], 'type')
                    if m_type is None:
                        m_type = getAttr(value[0], 'role')
                    m_ref = getAttr(value[0], 'ref')
                    reaction['P' + str(index+1) + 'Type'] = m_type
                    reaction['P' + str(index+1) + 'Ref'] = m_ref
                    self.molecules[m_ref]['type'] = m_type
            except:
                pass

            # reaction type
            if len(reactant_list) == 1 and len(product_list) == 1:
                reaction['type'] = 'A -> B'
            elif len(reactant_list) == 1 and len(product_list) == 2:
                reaction['type'] = 'A -> B + C'
            elif len(reactant_list) == 2 and len(product_list) == 1:
                reaction['type'] = 'A + B -> C'
            else:
                reaction['type'] = 'A + B -> C + D'

            # read MCRCMethod
            try:
                mcrc_method = getNode(node, 'MCRCMethod')
                # RRKM
                if getAttr(mcrc_method, 'name', from_text=True) == 'SimpleRRKM':
                    reaction['MCRCMethod'] = 'SimpleRRKM'
                    transition_state = getNode(node, 'transitionState')
                    m_type = getAttr(transition_state[0], 'type')
                    m_ref = getAttr(transition_state[0], 'ref')
                    reaction['TRef'] = m_ref
                    reaction['TType'] = m_type
                    if m_type is None:
                        m_type = getAttr(transition_state[0], 'role')
                    self.molecules[m_ref]['type'] = m_type
                # ILT parameters
                elif getAttr(mcrc_method, 'type') == 'MesmerILT' or getAttr(mcrc_method, 'name') == 'MesmerILT':
                    reaction['MCRCMethod'] = 'MesmerILT'
                    pre_exponential = getNode(node, 'preExponential')
                    reaction['preExponential'] = pre_exponential.text
                    reaction['preExponential_unit'] = getAttr(pre_exponential, 'units')

                    activation_energy = getNode(node, 'activationEnergy')
                    reaction['activationEnergy'] = activation_energy.text
                    reaction['activationEnergy_unit'] = getAttr(activation_energy, 'units')

                    n_infinity = getNode(node, 'nInfinity')
                    reaction['nInfinity'] = n_infinity.text
            except:
                pass
            reaction_list.append(reaction)

        return reaction_list


    def to_json(self):
        json_data = {}
        for el in self.root:
            if isinstance(el, etree._Comment):
                continue
            tag = etree.QName(el).localname
            if tag == 'title':
                json_data['title'] = el.text
            elif tag == 'moleculeList':
                json_data["moleculeList"] = self.read_molecule_list(el)
            elif tag == 'reactionList':
                json_data['reactionList'] = self.read_reaction_list(el)
            elif tag == 'conditions':
                json_data['conditions'] = self.read_condition(el)
            elif tag == 'modelParameters':
                json_data['modelParameters'] = self.read_model_parameters(el)
            elif tag == 'control':
                json_data['control'] = self.read_control(el)
        return json_data

# xml = MesmerXML(open('H:\mesmer\phenoxy_2.xml', 'r').read())
# data = xml.to_json()
# print data

