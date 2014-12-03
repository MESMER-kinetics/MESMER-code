from django.shortcuts import render
from django.http import HttpResponse
from django.template.loader import get_template
from django.template import Context
from lxml import etree
from xml.etree import cElementTree as ET
from xml.dom import minidom

from MesmerWeb import settings
from mesmerxml.xmlclass import MesmerXML

import os, subprocess, json, time, logging


# index for default request
def index(request):
    return render(request, 'html/index.html')


# parse a xml to json
def parsexml(request):
    if request.method == 'POST':
        try:
            xml_string = request.FILES['openfile'].read()
            xml = MesmerXML(xml_string)
            mesmer_json = xml.to_json()
            response = HttpResponse(json.dumps(mesmer_json))
            # response['Content-Type'] = "application/json; charset=utf-8"
            return response
        except:
            return HttpResponse('{"error": "error"}')
    else:
        return HttpResponse('{"error": "error"}')


# save xml file
def toxml(request):
    data = json.loads(request.body)
    t = get_template('xml/mesmer.xml')
    xml = t.render(Context(data))
    filename = "%s_%s.xml" % (str(time.time()), data.get('title'))
    filepath = os.path.abspath(os.path.join(settings.TMP_DIR, filename))
    with open(filepath, 'w') as destination:
        destination.write(pretty_xml(xml))
    return HttpResponse(filename)

def pretty_xml(xml_string):
    pretty = ''
    lines = xml_string.splitlines()
    for line in lines:
        line = line.rstrip()
        if line:
            pretty = pretty + line + '\n'
    return pretty


# download file
def download_xml(request, filename):
    path = os.path.join(settings.TMP_DIR, filename)
    return download(path)

def download(path):
    response = HttpResponse(read_file(path))
    response['Content-Disposition'] = "attachment;"
    response['Content-Type'] = "text/xml; charset=utf-8"
    return response


def read_file(path, buf_size=50000):
    f = open(path, "r")
    while True:
        block = f.read(buf_size)
        if block:
            yield block
        else:
            break
    f.close()
    os.remove(path)

def calculate(request):
    '''
    @param request:
    @return:
    '''
    log = logging.getLogger('test1')
    data = json.loads(request.body)
    t = get_template('xml/mesmer.xml')
    xml = t.render(Context(data))
    filename = ("%s_%s" % (str(time.time()), data.get('title'))).replace(' ', '_')
    filepath  = os.path.abspath(os.path.join(settings.TMP_DIR, filename+".xml"))
    output_filepath = os.path.abspath(os.path.join(settings.TMP_DIR, filename+"_out.xml"))
    with open(filepath, 'w') as destination:
        destination.write(pretty_xml(xml))
    destination.close()

    cal_process = subprocess.Popen(['mesmer', filepath, "-o", output_filepath], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    log.error("cal_process end")
    i = 0
    while(cal_process.poll() is None and i < 20):
        time.sleep(1)
        i += 1
    log.error("sleep end")
    if (cal_process.poll() is None):
        cal_process.kill()
    os.remove(filepath)
    log.error("begin to return")
    if os.path.exists(output_filepath):
        return HttpResponse('{"data": "%s"}' % (filename+"_out.xml"))
    else:
        err = cal_process.stderr.read()
        return HttpResponse('{"error": "%s"}' % repr(err)[1:-1])

def process_kill(p):
    p.kill()

def get_library(request):
    with open(os.path.join(settings.BASE_DIR, 'mesmerxml/molecule_library.json'), 'r') as f:
        library = f.read()
    return HttpResponse(library)