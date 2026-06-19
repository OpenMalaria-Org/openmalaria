#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of the openmalaria.tools package.
# For copyright and licensing information about this package, see the
# NOTICE.txt and LICENSE.txt files in its top-level directory; they are
# available at https://github.com/vecnet/openmalaria.tools
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License (MPL), version 2.0.  If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import argparse, sys, os.path, re
from datetime import date
#try:
    #import lxml.etree as ET
#except ImportError:
import xml.etree.ElementTree as ET

# document modes:
DM_NONE = 0
DM_HEAD = 1
DM_PARA = 2
DM_LIST = 3
DM_CODE = 4

class DocElts:
    def __init__(self, schema_file, ver, split):
        self.schema_file = schema_file
        self.ver = ver
        self.docname = 'schema-' + ver
        self.subdocs = { 'interventions' : self.docname + '-intervs', 'human' : self.docname + '-human' } if split else {}
        self.elements = []
        self.links = set()
    """Get the document name"""
    def name(self, elt):
        if elt.parent is None:
            return self.docname
        else:
            return self.subdocs.get(elt.name, elt.parent.docname)
    """Add an element"""
    def add(self, elt):
        self.elements.append(elt)
    """Create a new, unique, link"""
    def new_link(self, name):
        a = 'elt-' + name
        i = 1
        link = a
        while True:
            if not link in self.links:
                self.links.add(link)
                return link
            i += 1
            link = a + '-' + str(i)
    """Write the document"""
    def write(self, out_dir, commit):
        elt_dict = dict()
        for elt in self.elements:
            elts = elt_dict.setdefault(elt.docname, [])
            elts.append(elt)
        
        subdocs_rev = {v: k for k, v in self.subdocs.items()}
        
        for docname, elts in elt_dict.items():
            out_path = os.path.join(out_dir, docname + '.md')
            with open(out_path, 'w', encoding='UTF-8') as f_out:
                w = DocWriter(f_out, docname)
                filter = subdocs_rev.get(docname, None)
                w.header(self.schema_file, self.ver, filter, commit)
                for elt in elts:
                    elt.writedoc(w)
                w.finish()
        
        return self.docname

class DocWriter:
    """Writes formatted text to a file"""
    def __init__(self,f_out, docname=None):
        self.f = f_out
        self.mode = DM_NONE
        self.docname = docname
    
    """Switches mode, unless same and force=False. Adds a line-break if necessary."""
    def set_mode(self, mode, force):
        if not force and mode == self.mode:
            return
        # end previous mode:
        if self.mode == DM_CODE:
            self.line('```')
        # line break:
        if (mode == DM_HEAD and self.mode == DM_LIST) or \
            (mode == DM_PARA and self.mode in [DM_PARA, DM_LIST]) or \
            (mode == DM_LIST and self.mode == DM_LIST):
                self.line()
        # set new mode:
        self.mode = mode
    
    """wiki formatting functions: GitHub Markdown"""
    def heading(self, level, *args):
        self.set_mode(DM_HEAD, True)
        """level: 1 is biggest heading, 2 next, 5 smallest"""
        # Markdown allows underlining with = or -, but this is longer so don't do it
        #if level <= 2:
            #self.line(text)
            #char = '=' if level == 1 else '-'
            #self.line(char * max(len(text), 3))
        #else:
        self.line('#' * level, *args)
    
    """return an anchor"""
    def anchor(self, name):
        return '<a name="'+name+'"></a>'
    """return a link to target using text"""
    def link(self, docname, internname, text):
        """return text for a link, to be injected into a line"""
        if docname is None or docname == self.docname:
            assert internname is not None
            # intended for internal links in HTML documents, compatible with GitHub:
            target = '#' + internname
        elif internname is None:
            target = docname
        else:
            target = docname + '#' + internname
        return '['+text+']('+target+')'
    """return text wrapped with code to make it bold"""
    def bold(self, text):
        return '**' + text + '**'
    
    """Write a line, regardless of mode"""
    def line(self, *args):
        print(*args, file=self.f)
    """Write a line of paragraph text. Merges with previous line"""
    def p(self, *args):
        self.set_mode(DM_PARA, False)
        print(*args, file=self.f)
    """Write a line of paragraph text, with an explicit break. If no args, creates a new paragraph but no extra line."""
    def pn(self, *args):
        self.set_mode(DM_PARA, True)
        if len(args) > 0:
            print(*args, file=self.f)
    """write a line of a bulleted list, assuming new if last line wasn't bulleted"""
    def bulleted(self, *args):
        self.set_mode(DM_LIST, False)
        print('* ', *args, file=self.f)
    """Start a code block, with highlighting for lang"""
    def startcode(self, lang=''):
        self.set_mode(DM_CODE, True)
        print('```'+lang, file=self.f)
    """Write a line of code. Auto-switches if different from previous mode."""
    def code(self, *args):
        if DM_CODE != self.mode:
            # set new mode:
            self.line('```')
            self.mode = DM_CODE
        print(*args, file=self.f)
    
    """write header"""
    def header(self, schema_file, ver, filter, commit):
        self.heading(1, 'Schema', ver, 'documentation', ('' if filter is None else '('+filter+')'))
        self.p('Generated from: `'+schema_file+'`.')
        
        if commit:
            dt = date.today()
            datestr='{0}-{1}-{2}'.format(dt.year,dt.month,dt.day)
            self.pn('This state represents a prerelease for schema', schema_file, 'generated at '+datestr)
            self.p('of commit ['+commit+'](https://github.com/OpenMalaria-Org/openmalaria/commit/'+commit+') in branch _[develop](https://github.com/OpenMalaria-Org/openmalaria/tree/develop)_')
        
        self.pn('Key:')
        self.code('  abc           required (one)')
        self.code('[ def ]         optional (zero or one)')
        self.code('( ghi )*        any number (zero or more)')
        self.code('( jkl )+        at least one')
        self.code('( mno ){2,inf}  two or more occurrences')
        #self.p('<wiki:toc max_depth="3"/>')
    
    """Write appinfo"""
    def appinfo(self, appinfo):
        if 'units' in appinfo:
            self.p(self.bold('Units:'), appinfo['units'])
        if 'min' in appinfo:
            self.p(self.bold('Min:'), appinfo['min'])
        if 'max' in appinfo:
            self.p(self.bold('Max:'), appinfo['max'])
        #TODO: maybe also print any other items, excluding name? There shouldn't be any...
        self.pn()
    
    """End the document"""
    def finish(self):
        # end previous mode
        self.set_mode(DM_NONE, False)
        # write footer
        pass

def parse_appinfo(text):
    result = {}
    for pair in text.split(';'):
        parts = pair.split(':',1)
        if len(parts) == 2:
            result[parts[0].strip()] = parts[1].strip()
        else:
            assert len(pair.strip()) == 0 or die('bad appinfo:', text)
    return result

freq_symbs={
    (1,1) : (' ', ''),
    (0,1) : ('[', ']'),
    (0,'inf') : ('(', ')*'),
    (1,'inf') : ('(', ')+')
}
xsdns = 'http://www.w3.org/2001/XMLSchema'
xsdpre = '{' + xsdns + '}'
def replace_pre(name):
    if name is None:
        return None
    parts = name.split(':', 2)
    if len(parts) == 2 and parts[0] in nsmap:
        return '{' + nsmap[parts[0]] + '}' + parts[1]
    else:
        return name

class XSDType:
    """for built-in types (int, string, etc.)"""
    def __init__(self, name):
        self.name = name
        self.appinfo = {}
        self.doc = None
    """
    doc: document to add to
    stypes: schema types
    parent: element from which this element is visited
    """
    def collect_elements(self, doc, stypes, parent):
        pass
    def type_spec(self):
        return self.name
    def has_elts(self):
        return 1        # we don't have child elts but we do have a type
    def has_attrs(self):
        return False
    def write_attr_spec(self, w):
        pass
    def write_attr_doc(self, w):
        pass
    def write_elt_spec(self, w):
        w.code('    ' + self.name)
    def write_doc(self, w, docname):
        pass

class Node:
    """Base schema node"""
    def __init__(self, node):
        self.doc = None
        self.appinfo = {}
        #print('Parsing', node.tag, node.get('name'))
        for child in node:
            self.read_child(child)
    
    def read_child(self, child):
        if child.tag == xsdpre + 'annotation':
            doc_node = child.find(xsdpre + 'documentation')
            if doc_node is not None:
                self.doc = doc_node.text
            appinfo = child.find(xsdpre + 'appinfo')
            if appinfo is not None and appinfo.text is not None:
                self.appinfo = parse_appinfo(appinfo.text)
        else:
            die('unexpected child tag:', child.tag)

class ComplexType(Node):
    """Represent a type in the schema"""
    def __init__(self, node):
        self.attrs = []
        self.children = None # tree of child elements
        self.elements = [] # flat list of child elements
        self.base_name = None
        Node.__init__(self, node)
    
    def read_child(self, child):
        if child.tag == xsdpre + 'all':
            assert self.children is None
            self.children = ('all', self.read_elts(child))
        elif child.tag == xsdpre + 'sequence':
            assert self.children is None
            self.children = ('seq', self.read_elts(child))
        elif child.tag == xsdpre + 'choice':
            assert self.children is None
            self.children = ('choice', self.read_elts(child))
        elif child.tag == xsdpre + 'complexContent':
            self.read_content(child)
        elif child.tag == xsdpre + 'simpleContent':
            assert self.children is None
            self.children = 'simple'    # this is only really to say we can't have another element
            self.read_content(child)
        elif child.tag == xsdpre + 'attribute':
            self.attrs.append(Attribute(child))
        else:
            Node.read_child(self, child)
    
    def read_elts(self, parent):
        result = []
        for child in parent:
            if child.tag == xsdpre + 'element':
                elt = Element(child)
                result.append(elt)
                self.elements.append(elt)
            elif child.tag == xsdpre + 'sequence':
                result.append(('seq', self.read_elts(child)))
            elif child.tag == xsdpre + 'choice':
                result.append(('choice', self.read_elts(child)))
            elif child.tag == xsdpre + 'sequence':
                result.append(('sequence', self.read_elts(child)))
            else:
                die('unexpected child of <%s>: <%s>' %(parent.tag, child.tag))
        return result
    
    def read_content(self, parent):
        """Extension types"""
        for child in parent:
            if child.tag == xsdpre + 'extension':
                #Our extensions add attributes and sometimes a (seq) element
                # (where the base type is also seq or has no elements)
                assert self.base_name is None
                self.base_name = replace_pre(child.get('base'))
                for child2 in child:
                    self.read_child(child2)
    
    """
    doc: document to add to
    stypes: schema types
    parent: element from which this element is visited
    """
    def collect_elements(self, doc, stypes, parent):
        if hasattr(self, 'base_name'):
            for attr in self.attrs:
                attr.set_type(stypes)
            if self.base_name is None:
                self.base_type = None
            else:
                self.base_type = stypes.get(self.base_name)
                assert self.base_type is not None or die('unknown type',self.base_name)
                self.base_type.collect_elements(doc, stypes, parent)
            del self.base_name
        for child in self.elements:
            child.collect_elements(doc, stypes, parent)
    
    def has_attrs(self):
        return len(self.attrs) > 0
    def has_elts(self):
        have_elts = 2 if len(self.elements) > 0 else 0
        if have_elts == 0 and self.base_type is not None:
            return self.base_type.has_elts()
        return have_elts
    def write_elt_links(self, w):
        if self.base_type is not None:
            self.base_type.write_elt_links(w)
        for elt in self.elements:
            w.bulleted(w.link(elt.docname, elt.linkname, elt.name))
    def write_doc(self, w, docname):
        if self.doc is not None:
            w.heading(4, 'Documentation ('+docname+')')
            w.appinfo(self.appinfo)
            lines = self.doc.split('\n')
            for line in lines:
                w.line(line.lstrip())
        if self.base_type is not None:
            self.base_type.write_doc(w, 'base type')
    def write_attr_spec(self, w):
        if self.base_type is not None:
            self.base_type.write_attr_spec(w)
        for attr in self.attrs:
            if attr.use == 'optional' or attr.default:
                w.code('  [', attr.type_spec(),('] DEFAULT VALUE '+str(attr.default) if attr.default is not None else ']'))
            else:
                w.code('   ', attr.type_spec())
    def write_attr_doc(self, w):
        if self.base_type is not None:
            self.base_type.write_attr_doc(w)
        for attr in self.attrs:
            attr.writedoc(w)
    def write_elt_spec(self, w):
        if self.base_type is not None:
            self.base_type.write_elt_spec(w)
        elif self.children is not None and self.children != 'simple':
            self.write_elt_rec(self.children, w, 0)
    def write_elt_rec(self, node, w, depth):
        if isinstance(node, Element):
            symbs = freq_symbs.get(node.occurs, None)
            if symbs is None:
                symbs = ('(', '){'+str(node.occurs[0])+','+str(node.occurs[1])+'}')
            w.code('| '*depth + symbs[0], '<' + node.name, '... />', symbs[1])
        else:
            assert len(node) == 2
            mode = node[0]
            if mode == 'all':
                w.code('| '*depth + 'IN ANY ORDER:')
            elif mode == 'seq' or mode == 'sequence':
                w.code('| '*depth + 'IN THIS ORDER:')
            elif mode == 'choice':
                w.code('| '*depth + 'EXACTLY ONE OF:')
            else:
                assert False
            for child in node[1]:
                self.write_elt_rec(child, w, depth+1)

class Attribute(Node):
    def __init__(self, node):
        self.name = node.get('name')
        self.mode = None
        self.default = node.get('default')
        self.use = node.get('use') # optional or required
        if self.use is not None and self.use not in ('optional', 'required'):
            die('bad use specification on attribute:', node.attrib)
        if (self.default is None) == (self.use is None):
            if self.use is None:
                print('Warning: attribute',node.get('name'),\
                    'does not specify use or default; assuming optional', file=sys.stderr)
                self.use = 'optional'
            elif self.use == 'optional':
                print('Note: attribute', node.get('name'),\
                    'does not need to specify use="optional" in addition to default="..."', file=sys.stderr)
            else:
                print('Warning: attribute', node.get('name'),\
                    'specifies default="..." even though use="required"', file=sys.stderr)
        self.type_name = replace_pre(node.get('type')) # None if no 'type' attr
        Node.__init__(self, node)
    
    def read_child(self, child):
        if child.tag == xsdpre + 'simpleType':
            for child2 in child:
                assert child2.tag == xsdpre + 'restriction'
                # ignore 'base' attr (don't need)
                assert self.mode is None
                self.mode = 'enum'
                self.vals = []
                for child3 in child2:
                    assert child3.tag == xsdpre + 'enumeration'
                    assert len(child3) == 0
                    self.vals.append(child3.get('value'))
        else:
            Node.read_child(self, child)
    
    def set_type(self, stypes):
        if hasattr(self,'type_name'):
            if self.type_name is None:
                self.t = None
            else:
                self.t = stypes.get(self.type_name)
                assert self.t is not None or die('unknown type',self.type_name)
            del self.type_name
    
    def type_spec(self):
        rst=self.name + '='
        if self.mode == 'enum':
            sep='('
            for val in self.vals:
                rst += sep + '"' + val + '"'
                sep = ' or '
            rst += ')'
        else:
            rst += self.t.type_spec()
        return rst
    
    def writedoc(self, w):
        name = self.appinfo.get('name', self.name)
        w.heading(4, name)
        w.startcode('xml')
        w.code(self.type_spec())
        w.appinfo(self.appinfo)
        if self.default:
            w.p(w.bold('Default value:'), self.default)
        
        if self.doc is not None:
            w.pn()
            lines = self.doc.split('\n')
            for line in lines:
                line = line.lstrip();
                if len(line) > 0:
                    w.p(line.lstrip())

class Element(Node):
    """Represent an element type in the schema"""
    def __init__(self, node):
        self.name = node.get('name')
        self.elt_type = None
        self.type_name = replace_pre(node.get('type')) # None if no 'type' attr
        maxO = node.get('maxOccurs', 1)
        maxO = 'inf' if maxO == 'unbounded' else int(maxO)
        self.occurs = (int(node.get('minOccurs', 1)), maxO)
        Node.__init__(self, node)
    
    def read_child(self, child):
        if child.tag == xsdpre + 'complexType':
            assert self.elt_type is None and self.type_name is None or die('redundant type specification')
            self.elt_type = ComplexType(child)
        else:
            Node.read_child(self, child)
    
    def depth(self):
        return 0 if self.parent is None else 1 + self.parent.depth()
    
    """
    doc: document to add to
    stypes: schema types
    parent: element from which this element is visited
    """
    def collect_elements(self, doc, stypes, parent):
        if self in doc.elements:
            # Already in output document; use shortest bread-crumb trail
            if parent.depth() < self.parent.depth():
                self.parent = parent
                self.docname = doc.name(self)
        else:
            self.parent = parent
            self.docname = doc.name(self)
            self.linkname = doc.new_link(self.name)
            doc.add(self)
            
            if self.elt_type is None:
                self.elt_type = stypes.get(self.type_name)
            assert self.elt_type is not None or die('type not found:', self.type_name)
            self.elt_type.collect_elements(doc, stypes, self)
    
    def writedoc(self, w):
        w.heading(2, w.anchor(self.linkname), self.appinfo.get('name', self.name))
        w.p(self.breadcrumb(w, None))
        
        #w.heading(5, 'specification')
        w.startcode('xml')
        have_attrs = self.elt_type.has_attrs()
        have_elts = self.elt_type.has_elts()
        w.code('<'+self.name+('' if have_attrs else ('>' if have_elts>0 else '/>')))
        if have_attrs:
            self.elt_type.write_attr_spec(w)
            w.code('  >' if have_elts>0 else '  />')
        if have_elts>0:
            self.elt_type.write_elt_spec(w)
            w.code('</'+self.name+'>')
        
        #w.heading(3, 'Elements')
        if have_elts>1:
            self.elt_type.write_elt_links(w)
        
        if self.doc is not None:
            w.heading(4, 'Documentation (element)')
            w.appinfo(self.appinfo)
            lines = self.doc.split('\n')
            for line in lines:
                w.line(line.lstrip())
        self.elt_type.write_doc(w, 'type')
        
        if have_attrs:
            w.heading(3, 'Attributes')
            self.elt_type.write_attr_doc(w)
    
    def breadcrumb(self, w, parent):
        parent = self.parent if parent is None else parent
        r = '→ ' if parent is None else parent.breadcrumb(w, None) + ' → '
        return r + w.link(self.docname, self.linkname, self.name)

class FixedAttribute:
    def __init__(self, string):
        self.string = string
        self.use = 'required'
        self.default = None
    def set_type(self, stypes):
        pass
    def type_spec(self):
        return self.string
    def writedoc(self, w):
        pass

def translate(in_path, out_dir, schema_name, ver, split, commit=None):
    global nsmap
    
    # Read document:
    with open(in_path, 'r') as f_in:
        tree = ET.parse(f_in)
    root = tree.getroot()
    assert root.tag == xsdpre + 'schema'
    
    # Get namespace:
    targetns = root.get('targetNamespace')
    ompre = '' if targetns is None else '{' + targetns + '}'
    try:
        nsmap = root.nsmap
    except AttributeError:
        # default ElementTree implementation doesn't let us access this, so guess:
        nsmap = { 'xs' : xsdns }
        if targetns is not None:
            nsmap['om'] = targetns
    
    # Make a dictionary of types, built-in then custom:
    stypes = {
        xsdpre + 'boolean' : XSDType('boolean'),
        xsdpre + 'int' : XSDType('int'),
        xsdpre + 'integer' : XSDType('integer'),
        xsdpre + 'decimal' : XSDType('decimal'),
        xsdpre + 'double' : XSDType('double'),
        xsdpre + 'string' : XSDType('string')
    }
    # Also find the OM root element:
    omroot = None
    for child in root:
        if child.tag == xsdpre + 'complexType':
            stypes[ompre + child.get('name')] = ComplexType(child)
        elif child.tag == xsdpre + 'element':
            assert omroot is None or die('schema contains multiple root elements')
            omroot = Element(child)
        else:
            die('unexpected tag in schema:', child.tag)
    
    # handle namespace crap
    omroot.elt_type.attrs.append(FixedAttribute(
        'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'))
    if targetns is None: # pre 32
        omroot.elt_type.attrs.append(FixedAttribute(
            'xsi:noNamespaceSchemaLocation="'+schema_name+'"'))
    else:
        omroot.elt_type.attrs.append(FixedAttribute(
            'xmlns:om="http://openmalaria.org/schema/scenario_' + ver + '"'))
        omroot.elt_type.attrs.append(FixedAttribute(
            'xsi:schemaLocation="http://openmalaria.org/schema/scenario_'+ver+' '+schema_name+'"'))
    
    # Collect all linked elements into a list, by order visited:
    doc = DocElts(schema_name, ver, split)
    omroot.collect_elements(doc, stypes, None)
    
    # Write tthe output:
    return doc.write(out_dir, commit)

def die(*args):
    print('Error: ', *args, file=sys.stderr)
    sys.exit(1)

def make_key_list(v):
    try:
        return int(v[0])
    except:
        return 1000000

def main():
    # parse command-line arguments
    parser = argparse.ArgumentParser(description="""This tool converts one or more
        OpenMalaria schema files to a more readable documentation format.""")
    parser.add_argument('schema', metavar='SCHEMA', nargs='+',
        help="Schema file to be translated")
    parser.add_argument('-s','--split', action='store_true',
        help='Split output into multiple files. The exact split is hard-coded.')
    parser.add_argument('-i','--index', action='store_true',
        help='Generate an index')
    parser.add_argument('--mdbook', action='store_true', help='Generate SUMMARY.md for mdbook')
    parser.add_argument('-O','--out-dir', metavar='OUTDIR', action='store',
        help='Directory to output to. If not given the current directory is used.')
    parser.add_argument('-d','--develop', metavar='COMMIT' ,action='store',
        help='Add a comment that this is a pre-release based on commit COMMIT.')
    args = parser.parse_args()
    
    out_dir = '.' if args.out_dir is None else args.out_dir
    if not os.path.isdir(out_dir):
        die('Output directory doesn\'t exist (or isn\'t a directory):',out_dir)
    # why __rename not rename?
    
    # Translate files
    generated=[]
    for in_path in args.schema:
        in_name = os.path.basename(in_path)
        m = re.match('scenario_([0-9a-zA-Z_]+).xsd', in_name)
        assert m is not None or die('Expected schema files to have name scenario_*.xsd')
        ver = str(m.group(1)).replace('_', '-')
        
        print('Translating', in_path, file=sys.stderr)
        out_name = translate(in_path, out_dir, in_name, ver, args.split, commit=args.develop)
        major = ver.split('-')[0]
        generated.append((major, out_name, in_name))
    
    if args.index:
        path = os.path.join(out_dir, 'schema-index.md')
        print('Writing',path)
        with open(path, 'w') as f_out:
            w = DocWriter(f_out)
            w.heading(1, 'Generated Schema Documentation')
            w.p('This documentation was automatically generated from OpenMalaria schema (XSD) files.')
            w.p('Comments are welcome but be warned that edits will most likely be lost.')
            w.pn('')
            w.p('The script used to do this is')
            w.p('[generateDoc.py](https://github.com/OpenMalaria-Org/openmalaria/blob/main/util/generateDoc.py),')
            w.p('The command used to generate this page (after shell-expansion of wild-cards like `*`) is:')
            w.startcode('sh')
            w.code(' '.join(sys.argv))
            w.heading(2, 'Index')
            #TODO: instead of linking generated doc, we should search for everything matching schema-*.md
            for major, link, schema in sorted(generated, key = make_key_list, reverse=True):
                w.bulleted(w.link(link, None, 'Documentation for '+schema))
            w.finish()
    
    if args.mdbook:
        path=os.path.join(out_dir, 'intro.md')
        print('Writing',path)
        with open(path, 'w') as f_out:
            w = DocWriter(f_out)
            w.heading(1, 'Schema Documentation')
            w.p('This site displays OpenMalaria schema documentation in an easy-to-read form.')
            w.p('Links:')
            w.bulleted(w.link('schema-latest.md', None, 'Latest schema version'))
            w.bulleted(w.link('https://github.com/OpenMalaria-Org/openmalaria/wiki', None, 'OpenMalaria wiki'))
            w.bulleted(w.link('https://github.com/OpenMalaria-Org/openmalaria', None, 'OpenMalaria repository'))
            w.pn('')
            w.p('This documentation was automatically generated from OpenMalaria schema (XSD) files.')
            w.p('It is automatically kept up to date with schema files in the master branch of the repository.')
            dt = date.today()
            w.p('Contents were last updated on {0}-{1}-{2}.'.format(dt.year,dt.month,dt.day))
            w.finish()
        
        path=os.path.join(out_dir, 'SUMMARY.md')
        print('Writing',path)
        with open(path, 'w') as f_out:
            w = DocWriter(f_out)
            w.line(w.link('intro.md', None, 'Introduction'))
            w.pn('')
            w.bulleted(w.link('schema-latest.md', None, 'Latest version'))
            for major, link, schema in sorted(generated, key = make_key_list, reverse=True):
                w.bulleted(w.link(link+'.md', None,
                                  'Development version' if major == 'current' else 'Schema {}'.format(major)))
            w.finish()

if __name__ == "__main__":
    main()
