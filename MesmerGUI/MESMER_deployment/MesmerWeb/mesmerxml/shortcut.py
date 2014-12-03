__author__ = 'ciaosu'

xmlns = [
    'http://www.chem.leeds.ac.uk/mesmer',           # xmlns:me
    'http://www.xml-cml.org/schema',                # xmlns(default)
    'http://www.w3.org/2001/XMLSchema-instance',    # xmlns:cml
]


def getAttr(etree_node, attr, from_text=False):
    value = etree_node.get(attr)
    if not value:
        for ns in xmlns:
            value = etree_node.get("{%s}%s" % (ns, attr))
            if value:
                break
    if not value and from_text:
        value = etree_node.text
    try:
        value = value[value.index(":")+1:]
    except:
        pass
    return value


def getNodes(etree_node, tag):
    nodes = etree_node.xpath('.//%s' % tag)
    if not nodes:
        for ns in xmlns:
            nodes = etree_node.xpath('.//ns:%s' % tag, namespaces={'ns': ns})
            if nodes:
                break
    return nodes


def getNode(etree_node, tag):
    nodes = getNodes(etree_node, tag)
    try:
        return nodes[0]
    except:
        return None

def separate(str, ch):
    re = ''
    lines = str.splitlines()
    for line in lines:
        line = line.strip()
        if line:
            re = re + line + ch
    return re.rstrip()